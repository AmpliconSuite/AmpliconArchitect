"""Test the production fetch method without loading AA's reference repository.

Only the method is compiled from the source AST; tiny BAM/reference stand-ins
isolate read selection from AA initialization and its scientific dependencies.
"""
import ast
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace
import unittest
from zlib import crc32


def load_fetch():
    source = Path(__file__).resolve().parents[1] / "src" / "bam_to_breakpoint.py"
    tree = ast.parse(source.read_text())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef)
               and n.name == "bam_to_breakpoint")
    method = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == "fetch")
    namespace = {"crc32": crc32,
                 "hg": SimpleNamespace(chrLen={1: 1000}, chrNum=lambda c: 1)}
    exec(compile(ast.Module(body=[method], type_ignores=[]), str(source), "exec"), namespace)
    return namespace["fetch"]


class FakeBam:
    def __init__(self, records):
        self.records = records

    def fetch(self, chrom, start, stop):
        return (r for r in self.records if start <= r.position < stop)


class Reader:
    fetch = load_fetch()

    def __init__(self, records, ratio):
        self.bamfile = FakeBam(records)
        self.downsample_ratio = ratio


def paired_records():
    return [SimpleNamespace(query_name="read:" + str(i), position=i + offset + 1)
            for offset in (0, 200) for i in range(200)]


def probe():
    reads = Reader(paired_records(), 0.10083894600793622).fetch("chr1", 1, 400)
    selected = [(r.query_name, r.position) for r in reads]
    return {"count": len(selected),
            "sha256": hashlib.sha256(json.dumps(selected).encode()).hexdigest()}


class ReadDownsamplingTests(unittest.TestCase):
    def test_selection_is_independent_of_process_hash_seed(self):
        results = []
        for seed in (None, None, "0", "1", "2"):
            env = os.environ.copy()
            env.pop("PYTHONHASHSEED", None)
            if seed is not None:
                env["PYTHONHASHSEED"] = seed
            result = json.loads(subprocess.check_output(
                [sys.executable, str(Path(__file__).resolve()), "--probe"], env=env, text=True))
            self.assertGreater(result["count"], 0)
            self.assertLess(result["count"], 400)
            results.append(result)
        self.assertTrue(all(r == results[0] for r in results))

    def test_mates_and_overlapping_fetches_agree(self):
        reader = Reader(paired_records(), 0.10083894600793622)
        first = list(reader.fetch("chr1", 1, 200))
        second = list(reader.fetch("chr1", 201, 400))
        self.assertEqual([r.query_name for r in first], [r.query_name for r in second])
        overlap = list(reader.fetch("chr1", 100, 300))
        whole = list(reader.fetch("chr1", 1, 400))
        self.assertEqual(overlap, [r for r in whole if 100 <= r.position <= 300])

    def test_known_crc_bucket_and_existing_threshold(self):
        # IEEE CRC32("123456789") = 0xcbf43926; bucket 62.
        record = SimpleNamespace(query_name="123456789", position=1)
        self.assertEqual(list(Reader([record], 0.62).fetch("chr1", 1, 1)), [])
        self.assertEqual(list(Reader([record], 0.6201).fetch("chr1", 1, 1)), [record])
        self.assertEqual(list(Reader([record], 0).fetch("chr1", 1, 1)), [])

    def test_full_coverage_does_not_access_query_names(self):
        # No query_name attribute: the full-coverage path must bypass hashing.
        records = [SimpleNamespace(position=i) for i in range(1, 4)]
        self.assertEqual(list(Reader(records, 1).fetch("chr1", 1, 3)), records)


if __name__ == "__main__":
    if sys.argv[1:] == ["--probe"]:
        print(json.dumps(probe()))
    else:
        unittest.main()
