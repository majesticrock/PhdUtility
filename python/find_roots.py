import numpy as np
from scipy.optimize import brentq, newton, minimize_scalar

def find_all_roots(g, a=-1.0, b=1.0,
                   n_grid=201,
                   sign_tol=1e-12,
                   cluster_tol=1e-8,
                   refine_tol=1e-12,
                   min_search_width=1e-6):
    """
    Find all x in [a,b] solving g(x) == 0.
    
    Parameters
    ----------
    f : callable
        Function f(x). Must be vectorized or accept numpy arrays.
    E : float
        The target value.
    a, b : floats
        Interval endpoints (default [-1,1]).
    n_grid : int
        Number of grid points for initial scan (201 by default).
    sign_tol : float
        Tolerance for deciding exact zeros on grid.
    cluster_tol : float
        Distance below which two roots are considered the same.
    refine_tol : float
        Tolerance passed to root refiners.
    min_search_width : float
        Minimum width to use when creating local bracket for minimization.
    """
    # g(x) = f(x) - E
    xs = np.linspace(a, b, n_grid)
    gs = np.asarray(g(xs))

    roots = []

    # 1) Brackets with sign changes
    for i in range(len(xs)-1):
        g1, g2 = gs[i], gs[i+1]
        x1, x2 = xs[i], xs[i+1]

        # exact zero on grid
        if abs(g1) <= sign_tol:
            roots.append(x1)
            # still check interval for crossing to the right
            continue

        if g1 * g2 < 0:  # sign change -> bracket
            try:
                r = brentq(lambda x: g(x), x1, x2, xtol=refine_tol, rtol=refine_tol, maxiter=200)
                roots.append(r)
            except Exception:
                # fallback: use bisection manually if brentq fails
                lo, hi = x1, x2
                for _ in range(60):
                    mid = 0.5*(lo+hi)
                    gm = g(mid)
                    if abs(gm) < refine_tol:
                        lo = hi = mid
                        break
                    if (g(lo)) * gm <= 0:
                        hi = mid
                    else:
                        lo = mid
                roots.append(0.5*(lo+hi))

    # 2) Detect potential tangential (touching) roots where no sign change:
    # Find local minima of |g| on the grid
    abs_g = np.abs(gs)
    # candidate indices where a point is a local minimum on grid (or near zero)
    cand_idx = []
    for i in range(1, len(xs)-1):
        if abs_g[i] <= abs_g[i-1] and abs_g[i] <= abs_g[i+1]:
            # candidate if small enough or isolated trough
            if abs_g[i] < 1e-3 * (1.0 + np.max(abs_g)) or abs_g[i] < 1e-6:
                cand_idx.append(i)

    # also include any grid points with very small |g|
    very_small = np.where(abs_g <= sign_tol)[0]
    for ii in very_small:
        if ii not in cand_idx:
            cand_idx.append(ii)

    # refine each candidate by minimizing |g| in a small neighborhood
    for idx in cand_idx:
        # build a bracket around xs[idx]
        left = max(a, xs[idx] - (xs[1]-xs[0])*5)
        right = min(b, xs[idx] + (xs[1]-xs[0])*5)
        if right - left < min_search_width:
            # enlarge a bit
            left = max(a, xs[idx] - min_search_width)
            right = min(b, xs[idx] + min_search_width)

        # minimize |g(x)| in [left, right]
        res = minimize_scalar(lambda x: abs(g(x)), bounds=(left, right), method='bounded',
                              options={'xatol': 1e-8})
        if not res.success:
            continue
        x_min = res.x
        val = abs(g(x_min))
        if val <= max(sign_tol, 1e-8):
            # refine with Newton (or secant) starting from x_min
            try:
                r = newton(lambda x: g(x), x0=x_min, tol=refine_tol, maxiter=200)
            except Exception:
                # fallback: small bisection around x_min
                lo = max(a, x_min - 1e-6)
                hi = min(b, x_min + 1e-6)
                for _ in range(60):
                    mid = 0.5*(lo+hi)
                    if abs(g(mid)) < refine_tol:
                        lo = hi = mid
                        break
                    if (g(lo)) * (g(mid)) <= 0:
                        hi = mid
                    else:
                        lo = mid
                r = 0.5*(lo+hi)
            # check valid and inside interval
            if a - 1e-12 <= r <= b + 1e-12 and abs(g(r)) <= 1e-6:
                roots.append(r)

    # 3) endpoint checks (sometimes root exactly at boundary)
    for x_end in (a, b):
        if abs(g(x_end)) <= max(sign_tol, 1e-12):
            roots.append(x_end)

    # 4) deduplicate roots (cluster by proximity)
    if not roots:
        return np.array([])

    roots = np.array(roots)
    roots.sort()
    unique_roots = [roots[0]]
    for r in roots[1:]:
        if abs(r - unique_roots[-1]) > cluster_tol:
            unique_roots.append(r)

    # final filtering: ensure residual small
    final = []
    for r in unique_roots:
        if a - 1e-10 <= r <= b + 1e-10 and abs(g(r)) <= 1e-6:
            final.append(r)
    return np.array(final)
