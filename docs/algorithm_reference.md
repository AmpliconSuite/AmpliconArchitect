# AmpliconArchitect: Algorithm and Control Flow Reference

Jens Luebeck

---

## Overview

AmpliconArchitect reconstructs focal amplicon structures from WGS data. Given seed genomic intervals (regions of copy-number gain) and a BAM file, it:

1. Computes genome-wide coverage statistics
2. Extends seed intervals to find all connected amplified regions
3. Detects structural variant breakpoints from discordant read pairs
4. Segments the amplified region into copy-number intervals via mean-shift
5. Constructs a breakpoint graph (vertices = genomic positions with strand, edges = genomic segments or SVs)
6. Assigns copy numbers to all edges via convex optimization (MOSEK)
7. Decomposes the weighted graph into cycles representing amplicon structures

All source code is in `src/`. Entry point: `AmpliconArchitect.py`. Core processing: `bam_to_breakpoint.py`. Graph structures: `breakpoint_graph.py` + `abstract_graph.py`. Optimization: `mosek_solver.py`.

---

## Phase 0: Initialization

**File:** `AmpliconArchitect.py` (lines ~150–310)

### 0.1 Reference and BAM setup

- Reference annotations loaded from `$AA_DATA_REPO/<ref>/` via `ref_util` (`hg`): gene positions, oncogene list, centromeres, segmental duplications, mappability tracks, chromosome lengths/ordering.
- BAM opened with `pysam.Samfile()`.
- Seed BED file loaded as `interval_list`. Each interval is clipped to chromosome boundaries.

### 0.2 `bam_to_breakpoint` construction

The `bam_to_breakpoint` object (`b2b`) is the main stateful processing class. It is constructed once and reused for all amplicons.

**Key parameters set at construction:**

| Parameter | Source | Default |
|---|---|---|
| `pair_support_min` | `--pair_support_min` | None (auto-scaled; floor = 2 at 10x) |
| `foldback_pair_support_min` | `--foldback_pair_support_min` | None → inherits `pair_support` |
| `num_sdevs` | `--insert_sdevs` | 3.0 |
| `downsample` | `--downsample` | 0 (= auto 10x) |
| `sensitivems` | `--sensitivems` | False |
| `window_size` | hardcoded | 10,000 bp |
| `ms_window_size` | hardcoded | 10,000 bp |
| `mapping_quality_cutoff` | hardcoded | 5 |
| `breakpoint_mapping_quality_cutoff` | hardcoded | 20 |

### 0.3 Coverage statistics — `median_coverage()`

**File:** `bam_to_breakpoint.py` (~line 259)

**Algorithm:** Random genome-wide sampling

1. Select 1,000 non-conserved 10 kbp windows spread across all chromosomes, weighted by chromosome length, excluding centromeres and high-mappability-flagged regions.
2. Compute coverage in each window using `pysam` fetch with the downsampling ratio applied.
3. Compute median, mean, and std dev at two window sizes (10 kbp coarse, 300 bp fine) simultaneously.
4. Filter outlier windows: keep only windows where coverage is within 5× the median and > 0.
5. Estimate read-pair statistics by scanning a sample of reads: `read_length`, `insert_size` (median of proper-pair TLEN), `insert_std`, `percent_proper`.

**Key derived values:**

```
max_insert   = insert_size + num_sdevs * insert_std
min_insert   = insert_size - num_sdevs * insert_std
pair_support = max(
    round(wc_300_avg / 10.0 * (insert_size - read_length) / 2 / read_length / percent_proper),
    pair_support_min
)
```

When `pair_support_min` is auto (not user-set) and downsampling is active, the floor scales proportionally with the effective downsampled coverage:
```
effective_floor = max(2, round(2 * downsampled_coverage / 10.0))
pair_support = max(formula_result, effective_floor)
```

**Downsampling:** If `downsample > 0`, all subsequent BAM fetches probabilistically drop reads using `hash(read_name) % 100 / 100.0 < downsample_ratio`. The coverage stats are also scaled by `downsample_ratio` (or `sqrt(downsample_ratio)` for std dev fields).

**Caching:** Results are written to `$AA_DATA_REPO/coverage.stats` (tab-separated, one line per BAM path) and reused on subsequent runs for the same BAM and parameters.

**Output:** `basic_stats` tuple: `(wc_10k_median, wc_10k_avg, wc_10k_std, wc_300_median, wc_300_avg, wc_300_std, read_length, insert_size, insert_std, min_insert, max_insert, pair_support, percent_proper, num_sdevs, filesize)`

---

## Phase 1: Interval Discovery (EXPLORE mode)

**File:** `AmpliconArchitect.py` (~line 333), `bam_to_breakpoint.py`

For each seed interval, `interval_hops(ird, rdlist=all_ilist)` is called. The results are merged across all seeds to produce `irdgroups` — one group of intervals per amplicon.

### 1.1 `interval_hops()`

**Algorithm:** Priority-queue-based greedy graph exploration

The function treats the genome as a graph where nodes are amplified regions and edges are discordant read pairs connecting them.

```
State:
  seen_list   = intervals already explored
  unseen_list = min-heap of (-edge_weight, interval) tuples (max-weight priority)
  clist       = all discovered intervals (merged)

Loop (while unseen_list non-empty and len(seen_list) < 10):
  1. Pop highest-weight candidate interval from unseen_list
  2. Extend it via interval_extend()
  3. Detect discordant edges leaving it via interval_neighbors()
  4. For each outward edge:
       - Skip if target already covered by seen/clist
       - Skip if target region < 20 kbp and edge count < 2
       - Skip if edge has fewer than 2 discordant pairs
       - Otherwise: extend target region, push to unseen_list
                    with priority = discordant read count
  5. Add interval to seen_list and clist
```

### 1.2 `interval_extend()`

**Algorithm:** Exponential binary-search extension

From the seed boundary, extend outward (left and right independently) in geometrically increasing steps (1×, 2×, 4×, 8× of 10 kbp windows) checking whether each candidate extension is amplified. Stops at:
- Chromosome boundaries
- A region that fails `interval_amplified()` at increasing sizes
- A conserved or segdup region

### 1.3 `interval_amplified()`

Runs a quick mean-shift (`meanshift_segmentation()`) on the candidate region with bandwidth `hb=2`. A region is considered amplified if > 1/5 of windows have coverage above the amplification threshold:

- Standard mode: `coverage > max(arm_median + 3 * arm_std, 2.5 * genome_median)`
- Sensitive mode (`--sensitivems`): lower thresholds for high copy-number-variable or viral samples

### 1.4 Grouping

After all seeds are explored, the discovered intervals from all seeds are merged (`merge_clusters()`). Seeds whose discovered intervals overlap are grouped into the same amplicon. Each group becomes one `irdgroup` (an `interval_list`), which is processed as one amplicon.

---

## Phase 2: Per-Amplicon Reconstruction

**File:** `AmpliconArchitect.py` (~line 394), `bam_to_breakpoint.py`

For each `irdgroup`, `interval_filter_vertices(ilist, amplicon_name, runmode)` is called. This is the main reconstruction function; it internally orchestrates all remaining phases and writes the output files.

---

## Phase 3: Discordant Edge Detection

**File:** `bam_to_breakpoint.py` — `interval_discordant_edges()`

### 3.1 Read pair selection

Reads are fetched from the BAM across all intervals in `ilist`. A read pair is selected as **discordant** if:
- Both reads are mapped
- Not supplementary
- `is_proper_pair == False` (not expected concordant orientation/distance)
- Mapping quality ≥ `mapping_quality_cutoff` (5)
- Not an everted artifact: same chromosome, outward-facing orientation, distance < `max_insert` — these are candidates for tandem duplications and are handled separately as potential foldbacks

### 3.2 Clustering — Union-Find

**Data structure:** Each read pair is represented by an absolute genomic position pair `(absPos_read1, absPos_mate)`.

1. Sort read pairs by each end position independently (two sorted lists).
2. For each read pair `v`, find all other pairs `u` where:
   - `|u.pos1 - v.pos1| < max_insert - read_length` AND
   - `|u.pos2 - v.pos2| < max_insert - read_length`
   This is done efficiently via binary search on the sorted lists, then intersecting.
3. All pairs in the intersection are connected in a Union-Find structure (union by rank).
4. Any cluster with ≥ `pair_support` members is retained.

**Breakpoint vertex computation:** For each cluster pair `(c1, c2)`:
- Forward-strand end: `pos = max(read.reference_end - 1)` across reads in cluster
- Reverse-strand end: `pos = min(read.reference_start)` across reads in cluster  
- Vertex: `breakpoint_vertex(chrom, pos, strand)`

### 3.3 Foldback detection

A candidate breakpoint is a **foldback** if:
- Both endpoints on same chromosome and same strand
- `|pos1 - pos2| ≤ read_length`

For foldbacks, reads are partitioned into *inverted* (both reads start at same position, or SA-tag overlaps primary) vs. *non-inverted*. Only the non-inverted read count is compared against `foldback_pair_support_min = max(foldback_pair_support_min, ps)`.

### 3.4 Edge refinement — `refine_discordant_edge()`

For each candidate breakpoint edge `(v1, v2)`:
1. Fetch all reads in `[v1.pos ± max_insert] × [v2.pos ± max_insert]`.
2. Remove multi-mapped reads (appear > 1 time in either region).
3. Compute **k-mer homology**: count 10-mers shared between 200 bp flanks on each side of the breakpoint (using set intersection). Stored in `edge.hom`.
4. Extract soft-clipped junction sequence from CIGAR ops. Stored in `edge.hom_seq`.

**Repeat filtering:** By default, edges are filtered if both endpoints fall in segmental duplications or low-mappability regions (configurable).

**Output:** `eilist` = list of `(breakpoint_edge, pair_count)` tuples

---

## Phase 4: Copy-Number Segmentation

**File:** `bam_to_breakpoint.py` — `meanshift_refined()`

### 4.1 Two-stage mean-shift segmentation

**Stage 1 — Coarse (10 kbp windows):**

`meanshift_segmentation()` is run on the amplicon interval with:
- Window size: 10,000 bp
- Significance threshold: `pvalue = 0.0027` (corresponds to ~3 sigma; Bonferroni-corrected for multiple comparisons)
- Kernel bandwidths tried in sequence: `hb ∈ {2, 5, 10, 50, 100}`

This produces coarse segment boundaries `segments0`.

**Stage 2 — Fine refinement (300 bp windows):**

For each transition in `segments0`, a ±3-window neighbourhood is re-segmented at 300 bp resolution with `pvalue = 0.05`. Segment boundaries are updated if:
- A fine boundary is within 10 kbp of the coarse boundary
- The copy-number change direction is consistent (`cn_diff_coarse * cn_diff_fine > 0`)
- The magnitude ratio is between 0.5 and 2.0

Boundaries passing these checks are marked `start_refined = True` / `end_refined = True`.

### 4.2 Mean-shift algorithm — `meanshift_segmentation()`

The implementation is a 1D adaptive-bandwidth mean-shift applied to a coverage vector.

**Coverage vector:** `c[i]` = mean coverage in window `i` (10 kbp or 300 bp)

**Adaptive kernel bandwidth:**
```
h(i) = sqrt(c[i] / rd_global) * h0    if c[i] >= rd_global / 4
      = h0 / 2                          if c[i] < rd_global / 4
```
where `h0` is the base bandwidth (selected from `{2, 5, 10, 50, 100}`) and `rd_global` is the global median coverage.

**Mean-shift vector** at window `wi`:
```
dfi[wi] = sum over wj in [-n, +n]:
    wj * exp(-0.5 * wj^2 / hb^2) * exp(-0.5 * (c[wi+wj] - c[wi])^2 / h(wi)^2)
```

**Segment detection:**
1. Find zero-crossings of `dfi` (sign changes from + to −).
2. Validate each boundary with a two-sample t-test between adjacent segments. Only boundaries with `p < pvalue` are retained.
3. Multi-scale refinement loop: run coarse → fine → coarse to stabilize boundaries.

**Output:** List of intervals with `info = {'cn': copy_number, 'start_refined': bool, 'end_refined': bool}`

---

## Phase 5: Breakpoint Graph Construction

**File:** `bam_to_breakpoint.py` — `interval_filter_vertices()`  
**Graph classes:** `breakpoint_graph.py`

### 5.1 Graph model

**Vertices** (`breakpoint_vertex`): A genomic position with strand — `(chrom, pos, strand)` where `strand ∈ {+1, −1}`. Vertices are stored in a hash map keyed by string representation `chrom:pos+/-`. One special **source vertex** `(chrom, -1, -1)` represents connections to the "rest of the genome" (reads with mates outside the amplicon).

**Edges** have four types:

| Type | Connects | Represents |
|---|---|---|
| `sequence` | `chrom:pos−` → `chrom:pos+` (same chrom, opposite strands) | A contiguous genomic segment |
| `concordant` | Adjacent positions, opposite strands (distance = 1 bp) | The reference genome connection between two adjacent segments |
| `discordant` | Any two vertices | A structural variant breakpoint |
| `source` | Source vertex → any vertex | Reads spanning the amplicon boundary |

A sequence edge from `v1 = chrom:start-` to `v2 = chrom:end+` represents the interval `[start, end]` on the chromosome. Traversing this edge in a cycle means including that genomic segment.

### 5.2 Graph construction

1. **Vertices from CN segments:** For each mean-shift segment boundary, create one `pos+` and one `pos−` vertex. Adjacent segments share a concordant edge at their boundary (breakpoint position ± 1 bp).

2. **Sequence edges:** Connect `start−` → `end+` for each segment. Weight (read count) = sum of read coverage across the segment.

3. **Source edges:** For each interval boundary (amplicon start/end), count proper read pairs with one mate inside and one outside. The outside-mate count becomes the source edge weight (`koe`).

4. **Discordant edges:** Each entry in `eilist` is added as a discordant edge between the two corresponding vertices. If one endpoint falls outside all intervals, it connects to the source vertex. Weight = pair count from clustering step.

5. **Concordant edges:** For each pair of adjacent segment boundaries (a `pos+` and the next `pos−`), count proper read pairs spanning both positions. Added if count ≥ `pair_support`.

**Output:** A `breakpoint_graph` object with weight dicts: `koe` (source edges), `kbpe` (breakpoint edges), `kce` (concordant edges), and coverage per sequence edge.

---

## Phase 6: Copy-Number Optimization

**File:** `mosek_solver.py`

### 6.1 Problem formulation

Let:
- `n` = number of sequence edges
- `m` = number of breakpoint + concordant + source edges
- `x ∈ ℝ^(n+m)` = unknown copy counts for each edge

**Objective (MOSEK v10+):**
```
minimize   c^T * x  −  Σ_i  f_i * log(x_i)
subject to  A * x = 0      (flow conservation)
            x ≥ 0
```

This is a maximum-likelihood formulation: the log terms reward solutions where copy count `x_i` is proportional to observed read count `f_i`. The linear term `c^T * x` penalises high copy numbers (acts as a regulariser, scaled by genomic length).

**Older MOSEK (v8):** Solves `minimize c^T*x − Σ f_i * log(g_i * x_i + h_i)` with normalising constants `g_i` and `h_i` that account for the expected number of reads per copy per edge.

### 6.2 Coefficient derivation

Let `C = wc_300_avg / 2` (estimate of haploid coverage).

For **sequence edge** `i` with length `l_i` and observed read count `k_i`:
```
f_i = k_i                          (observed read count)
g_i = C * l_i / read_length        (expected reads per copy)
c_i = g_i                          (prior penalty)
```

For **breakpoint/concordant/source edge** `j` with observed pair count `k_j`:
```
f_j = k_j
g_j = C * max_insert / 2 / read_length
c_j = g_j
```

### 6.3 Flow conservation constraints

The matrix `A` encodes Eulerian flow balance: for each sequence edge `i`, the copy count of that sequence edge must equal the sum of copy counts of all non-sequence edges incident to each of its two endpoints.

This gives `2n` equality constraints (one per sequence-edge endpoint):
```
x_seq[i] = Σ x_bp[j]  for all j incident to v1(i)
x_seq[i] = Σ x_bp[j]  for all j incident to v2(i)
```

### 6.4 Solver dispatch

| MOSEK version | API used | Notes |
|---|---|---|
| 8 | `scopt` (self-concordant optimizer) | Includes `g_i, h_i` normalisation |
| 9 | Fusion API, primal exponential cone | Dropped `h_i` |
| ≥ 10 | ACC (affine conic constraints) | `minimize c^T*x − Σ f_i*log(x_i)`, `x > 0` implicit |

Tolerance (v10+): `intpnt_co_tol_near_rel = 1e5`. On solver failure, falls back to returning `coeff_c` as copy numbers and logs the inputs to a JSON file.

**Output:** `res[0:n]` = sequence edge copy counts, `res[n:n+m]` = breakpoint/concordant/source edge copy counts.

---

## Phase 7: Cycle Decomposition

**File:** `breakpoint_graph.py` — `cycle_decomposition(w, s)`

### 7.1 Problem

Given the edge weight dict `w` (copy counts from MOSEK) and the source vertex `s`, decompose the graph into a set of cycles (or paths through the source vertex) that cover the amplicon content. Each cycle corresponds to one structural form of the amplicon.

**Stopping criterion:** Continue until ≥ 80% of total amplicon content (sum of `sequence_edge_length × copy_count`) is explained.  
**Minimum cycle weight:** 0.1 copies.

### 7.2 Thickest-cycle-first greedy algorithm

```
w2 = copy of w
cycle_list = []

while max(w2.values()) > 0.1:
    For each edge e in w2 (sorted by weight, highest first):
        tc, tcw = thickest_cycle(e, w2)  # find max-weight cycle through e
    
    Select tc = cycle with maximum tcw (minimum weight along cycle)
    For each edge in tc: w2[edge] -= tcw
    cycle_list.append((cycle_number, tcw, tc))
    cycle_number += 1
    if amplicon_content_covered >= 0.80 * total_amplicon_content:
        break
```

### 7.3 `thickest_cycle(e, w)` — max-weight cycle search

This is a modified Dijkstra's algorithm that finds the maximum-weight cycle containing edge `e`.

**Graph traversal rule:** In a valid cycle, breakpoint/concordant/source edges and sequence edges must alternate. The algorithm therefore always crosses a breakpoint edge then its adjacent sequence edge together (one "hop" = one non-sequence edge + one sequence edge).

**State:** `hdict[v] = (best_weight, path_edges, predecessor, seen_edge_set)`  
**Heap:** min-heap of `(-weight, vertex)` (negated for max-weight search)

```
Start: push (edge_weight[e], v1_of_e) onto heap

At each vertex v1 popped:
  If v1 == starting_vertex and already seen: cycle found, stop

  For each non-sequence edge e incident to v1:
    v2 = other endpoint of e
    
    If v2 == source:
      new_weight = min(current_weight, w[e])  [half if e already in path]
      Update hdict[source] if better
    
    Else:
      Find the unique sequence edge se adjacent to v2
      v3 = other endpoint of se
      new_weight = min(current_weight, w[e])  [halved if e or se repeated]
      Update hdict[v3] if better, recording [e, se] as path step
    
    Push (new_weight, v3) to heap
```

**Repeated edge penalty:** If a breakpoint edge or sequence edge is already in the current path's `seen_edge_set`, its contribution to the minimum weight is halved (allows reuse but penalises it, preventing infinite loops).

**Cycle extraction:** Backtrack through `hdict` from the start vertex to reconstruct the edge list.

### 7.4 Cycle orientation and output

After finding a cycle:
1. If the cycle passes through the source vertex, reorder to start after the source.
2. Ensure the first sequence edge traverses in ascending genomic position (flip entire cycle if not).
3. Rotate the cycle so it does not start with a concordant edge.

Each cycle is logged to `_cycles.txt` in the format:
```
Segment  {id}  {chrom}  {start}  {end}
...
Cycle={id};Copy_count={weight:.2f};Segments={seg_id}{+/-},...
```

Orientation: `+` = sequence edge traversed in forward direction (ascending position), `−` = reverse.

---

## Output Files

### `_graph.txt`

```
SequenceEdge: StartPosition EndPosition PredictedCopyCount AverageCoverage Size NumberOfReadsMapped
BreakpointEdge: StartPosition->EndPosition PredictedCopyCount NumberOfReadsMapped HomologySizeinGraph MicrohomologyOrInsertionSequence
```

Edge types in the file: `sequence`, `concordant`, `discordant`, `source`.

### `_cycles.txt`

```
Segment  {n}  {chrom}  {start}  {end}
...
Cycle={n};Copy_count={float};Segments={seg+/-,...}
```

Segment IDs are 1-indexed and ordered by genomic position within the amplicon. A segment followed by `−` means the reverse complement of that region is included in the cycle at that position.

### `_cnseg.txt` (intermediate, cached)

Tab-separated: `chrom  start  end  copy_number  start_refined  end_refined`

### `_edges.txt` (intermediate, cached)

One breakpoint edge per line: `edge_string  pair_count  homology_bp  junction_sequence`

---

## Key Thresholds and Parameters Summary

| Parameter | Default | Controlled by |
|---|---|---|
| `pair_support` | auto (≥2 at 10x) | Coverage formula + `--pair_support_min` |
| `foldback_pair_support_min` | = `pair_support` | `--foldback_pair_support_min` |
| `max_insert` | `insert_size + 3σ` | `--insert_sdevs` |
| Coarse window size | 10,000 bp | Hardcoded |
| Fine window size | 300 bp | Hardcoded |
| Mean-shift bandwidths | {2, 5, 10, 50, 100} windows | Hardcoded |
| Coarse segmentation p-value | 0.0027 | Hardcoded |
| Fine segmentation p-value | 0.05 | Hardcoded |
| Amplification threshold (std) | arm median + 3σ | `--sensitivems` changes this |
| Cycle stop threshold | 80% content covered | Hardcoded |
| Min cycle copy count | 0.1 | Hardcoded |
| k-mer homology k | 10-mer, ±100 bp span | Hardcoded |
| Mapping quality (general) | 5 | Hardcoded |
| Mapping quality (breakpoints) | 20 | Hardcoded |
| Max exploration hops per seed | 10 | Hardcoded (`interval_hops`): caps how many connected intervals can be discovered by following discordant edges outward from a single seed; does not limit the number of seeds in the BED file |
| Haploid coverage estimate | `wc_300_avg / 2` | Derived |
| MOSEK near-rel tolerance | 1e5 | Hardcoded |
