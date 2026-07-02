# Interface to Clarabel for AmpliconArchitect copy-number optimization.
import logging

import numpy as np
from scipy import sparse


clarabel_logger = logging.getLogger('Clarabel')


def _flow_residual(x, asub, aval):
    max_residual = 0.0
    for cols, vals in zip(asub, aval):
        residual = sum(x[col] * val for col, val in zip(cols, vals))
        max_residual = max(max_residual, abs(residual))
    return max_residual


def _objective(x, coeff_c, coeff_f):
    if np.any(x <= 0):
        return np.inf
    return float(np.dot(coeff_c, x) + np.dot(coeff_f, np.log(x)))


def _valid_solution(x, asub, aval, coeff_c, coeff_f):
    if not np.all(np.isfinite(x)):
        return False
    if np.any(x < -1e-7):
        return False
    if not np.isfinite(_objective(np.maximum(x, 1e-300), coeff_c, coeff_f)):
        return False
    scale = max(1.0, float(np.max(np.abs(x)))) if len(x) else 1.0
    if _flow_residual(x, asub, aval) > 1e-5 * scale:
        return False
    return True


def _settings(clarabel, relaxed=False, alternate_preprocessing=False):
    settings = clarabel.DefaultSettings()
    settings.verbose = False
    settings.max_threads = 1
    settings.max_iter = 1000

    if relaxed:
        settings.tol_feas = 1e-6
        settings.tol_gap_abs = 1e-6
        settings.tol_gap_rel = 1e-6
        settings.tol_infeas_abs = 1e-6
        settings.tol_infeas_rel = 1e-6
        settings.tol_ktratio = 1e-6

    if alternate_preprocessing:
        settings.presolve_enable = False
        settings.equilibrate_enable = False

    return settings


def call_clarabel(n, m, asub, aval, coeff_c, coeff_f):
    """Solve the MOSEK v10 AA copy-number model with Clarabel.

    The AA model is:

        minimize    c^T x + f^T log(x)
        subject to  A x = 0

    where ``coeff_f`` contains non-positive read-count coefficients. Clarabel
    solves conic problems in the form ``P/2*x^2 + q^T*x`` subject to
    ``A*x + s = b, s in K``. We introduce one auxiliary variable ``t_i`` for
    each edge and constrain ``(t_i, 1, x_i)`` to the exponential cone. The
    objective ``c^T x + f^T t`` then makes ``t_i = log(x_i)`` at optimum.
    """
    try:
        import clarabel
    except ImportError as e:
        raise ImportError("Clarabel solver requested but the 'clarabel' Python package is not installed") from e

    clarabel_logger.info("Beginning Clarabel call")

    num_edges = n + m
    num_vars = 2 * num_edges
    num_flow_cons = 2 * n
    num_exp_cons = 3 * num_edges

    q = np.array(list(coeff_c) + list(coeff_f), dtype=float)
    p = sparse.csc_matrix((num_vars, num_vars))
    b = np.zeros(num_flow_cons + num_exp_cons, dtype=float)

    rows = []
    cols = []
    data = []

    for row_idx, (row_cols, row_vals) in enumerate(zip(asub, aval)):
        for col, val in zip(row_cols, row_vals):
            rows.append(row_idx)
            cols.append(col)
            data.append(val)

    for i in range(num_edges):
        base = num_flow_cons + 3 * i

        rows.append(base)
        cols.append(num_edges + i)
        data.append(-1.0)

        b[base + 1] = 1.0

        rows.append(base + 2)
        cols.append(i)
        data.append(-1.0)

    a_mat = sparse.csc_matrix((data, (rows, cols)), shape=(num_flow_cons + num_exp_cons, num_vars))
    cones = [clarabel.ZeroConeT(num_flow_cons)]
    cones.extend(clarabel.ExponentialConeT() for _ in range(num_edges))

    attempts = [
        ("default", _settings(clarabel)),
        ("relaxed_tolerances", _settings(clarabel, relaxed=True)),
        ("relaxed_no_presolve_equilibration", _settings(clarabel, relaxed=True, alternate_preprocessing=True)),
    ]

    for attempt_name, settings in attempts:
        clarabel_logger.info("Clarabel attempt: {}".format(attempt_name))
        solver = clarabel.DefaultSolver(p, q, a_mat, b, cones, settings)
        solution = solver.solve()
        status = str(solution.status)
        x = np.array(solution.x[:num_edges], dtype=float)

        if status in ("Solved", "AlmostSolved") and _valid_solution(x, asub, aval, coeff_c, coeff_f):
            if attempt_name != "default":
                logging.warning("Clarabel succeeded with fallback attempt '{}', status {}.".format(attempt_name, status))
            return x

        clarabel_logger.warning("Clarabel attempt '{}' failed validation with status {}.".format(attempt_name, status))

    logging.error("Failed to solve to optimality with Clarabel after all retry attempts.")
    logging.error("Due to failure to solve, copy numbers reported by AA will not be reliable!")
    return np.array(coeff_c)
