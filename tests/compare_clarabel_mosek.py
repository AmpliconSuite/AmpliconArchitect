#!/usr/bin/env python
"""Compare AA's MOSEK solve against a Clarabel solve on the same inputs.

This is a development probe, not a committed runtime backend. It accepts the
JSON files written by ``mosek_solver.save_mosek_input`` and can also generate a
few small AA-shaped synthetic problems.
"""
import argparse
import json
import math
import os
import sys
import time

import numpy as np


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC_DIR = os.path.join(REPO_ROOT, "src")
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

import mosek_solver  # noqa: E402
import clarabel_solver  # noqa: E402


def objective(x, coeff_c, coeff_f):
    """Return the AA v10 objective value c'x + f't where t=log(x)."""
    x = np.asarray(x, dtype=float)
    coeff_c = np.asarray(coeff_c, dtype=float)
    coeff_f = np.asarray(coeff_f, dtype=float)
    if np.any(x <= 0):
        return math.inf
    return float(np.dot(coeff_c, x) + np.dot(coeff_f, np.log(x)))


def max_flow_residual(x, asub, aval):
    residuals = []
    for cols, vals in zip(asub, aval):
        residuals.append(sum(x[col] * val for col, val in zip(cols, vals)))
    return max(abs(r) for r in residuals) if residuals else 0.0


def load_problem(path):
    with open(path) as handle:
        data = json.load(handle)
    return (
        data["n"],
        data["m"],
        data["asub"],
        data["aval"],
        data["coeff_c"],
        data["coeff_f"],
        data.get("coeff_g"),
        data.get("const_h"),
    )


def synthetic_problems():
    problems = []

    # One sequence edge and one non-sequence edge with identical endpoint flow.
    problems.append(
        (
            "one_cycle",
            1,
            1,
            [[0, 1], [0, 1]],
            [[1.0, -1.0], [1.0, -1.0]],
            [20.0, 20.0],
            [-120.0, -80.0],
            [20.0, 20.0],
            [0.0001, 0.0001],
        )
    )

    # Two sequence edges joined by three non-sequence edges. The duplicate
    # endpoint constraints mirror the shape produced by bam_to_breakpoint.py.
    problems.append(
        (
            "two_segment_chain",
            2,
            3,
            [[0, 2], [0, 3], [1, 3], [1, 4]],
            [[1.0, -1.0], [1.0, -1.0], [1.0, -1.0], [1.0, -1.0]],
            [30.0, 28.0, 18.0, 16.0, 15.0],
            [-180.0, -165.0, -95.0, -90.0, -85.0],
            [30.0, 28.0, 18.0, 16.0, 15.0],
            [0.0001] * 5,
        )
    )

    # Two sequence edges with two non-sequence paths. This gives the optimizer
    # freedom to split flow between the non-sequence edges.
    problems.append(
        (
            "branched_pair",
            2,
            2,
            [[0, 2, 3], [0, 2, 3], [1, 2, 3], [1, 2, 3]],
            [[1.0, -1.0, -1.0]] * 4,
            [35.0, 33.0, 14.0, 21.0],
            [-210.0, -198.0, -55.0, -115.0],
            [35.0, 33.0, 14.0, 21.0],
            [0.0001] * 4,
        )
    )

    return problems


def compare_problem(name, problem, rtol, atol):
    n, m, asub, aval, coeff_c, coeff_f, coeff_g, const_h = problem

    start = time.time()
    mosek_res = np.asarray(
        mosek_solver.call_mosek(n, m, asub, aval, coeff_c, coeff_f, coeff_g, const_h),
        dtype=float,
    )
    mosek_seconds = time.time() - start

    start = time.time()
    clarabel_res = clarabel_solver.call_clarabel(n, m, asub, aval, coeff_c, coeff_f)
    clarabel_seconds = time.time() - start

    abs_diff = np.abs(mosek_res - clarabel_res)
    denom = np.maximum(np.abs(mosek_res), atol)
    rel_diff = abs_diff / denom
    ok = np.allclose(mosek_res, clarabel_res, rtol=rtol, atol=atol)

    print("CASE {}".format(name))
    print("  status: {}".format("PASS" if ok else "FAIL"))
    print("  mosek_seconds: {:.6f}".format(mosek_seconds))
    print("  clarabel_seconds: {:.6f}".format(clarabel_seconds))
    print("  comparison_rtol: {:.6g}".format(rtol))
    print("  comparison_atol: {:.6g}".format(atol))
    print("  max_abs_diff: {:.6g}".format(float(np.max(abs_diff))))
    print("  max_rel_diff: {:.6g}".format(float(np.max(rel_diff))))
    print("  mosek_objective: {:.12g}".format(objective(mosek_res, coeff_c, coeff_f)))
    print("  clarabel_objective: {:.12g}".format(objective(clarabel_res, coeff_c, coeff_f)))
    print("  mosek_flow_residual: {:.6g}".format(max_flow_residual(mosek_res, asub, aval)))
    print("  clarabel_flow_residual: {:.6g}".format(max_flow_residual(clarabel_res, asub, aval)))
    print("")

    return ok


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="*", help="mosekinput-*.json files to compare")
    parser.add_argument("--synthetic", action="store_true", help="run built-in synthetic cases")
    parser.add_argument("--rtol", type=float)
    parser.add_argument("--atol", type=float, default=1e-4)
    args = parser.parse_args()

    cases = []
    if args.synthetic or not args.inputs:
        for item in synthetic_problems():
            name, *problem = item
            cases.append((name, tuple(problem)))

    for path in args.inputs:
        cases.append((path, load_problem(path)))

    failed = 0
    rtol = args.rtol
    if rtol is None:
        rtol = 1e-4

    for name, problem in cases:
        if not compare_problem(name, problem, rtol, args.atol):
            failed += 1

    if failed:
        print("{} comparison case(s) failed".format(failed), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
