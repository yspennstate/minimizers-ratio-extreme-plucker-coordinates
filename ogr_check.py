"""Numerical checks of Ogranovich Conjecture 1.2 (k >= 3) and of the extreme-minor formula.
V0: span of Re/Im of (w^0..w^{n-1}) for the k eigenvalues w = exp(i pi (k-1-2r)/n), r = 0..k-1, of the signed shift."""
import itertools, math, sys
import numpy as np
from scipy.optimize import minimize

rng = np.random.default_rng(1)


def v0(k, n):
    rows = []
    ws = [np.exp(1j * math.pi * (k - 1 - 2 * r) / n) for r in range(k)]
    for r in range(k):
        v = ws[r] ** np.arange(n)
        if ws[r].imag > 1e-12:
            rows += [v.real, v.imag]
        elif abs(ws[r].imag) < 1e-12:
            rows.append(v.real)
    M = np.array(rows[:k])
    assert np.linalg.matrix_rank(M) == k
    return M


def minors(M, subs):
    return np.array([np.linalg.det(M[:, list(I)]) for I in subs])


def E(M, subs):
    d = minors(M, subs)
    if d.min() <= 0:
        return np.inf
    return d.max() / d.min()


out = []
for k, n in [(3, 5), (3, 6), (3, 7), (3, 8), (4, 7), (4, 8)]:
    subs = list(itertools.combinations(range(n), k))
    M0 = v0(k, n)
    d0 = minors(M0, subs)
    if d0.min() < 0:
        d0 = -d0
        M0[0] = -M0[0]
    P = np.array([np.prod([math.sin(math.pi * (b - a) / n) for a, b in itertools.combinations(I, 2)]) for I in subs])
    ratio_err = np.ptp(d0 / P) / np.mean(d0 / P)
    Bset = tuple(sorted(math.floor(j * n / k) for j in range(k)))
    e0 = E(M0, subs)
    formula = P[subs.index(Bset)] / P[subs.index(tuple(range(k)))]
    # local + random search for a smaller E; parametrize near V0 and from random positive starts
    best = np.inf
    used = 0
    f = lambda x: E(x.reshape(k, n), subs)
    for trial in range(25):
        if trial < 15:
            x0 = (M0 + 0.05 * (trial + 1) * rng.standard_normal(M0.shape)).ravel()
        else:  # random totally positive start: points on a moment curve with random nodes
            t = np.sort(rng.uniform(0, 1, n)) + np.arange(n)
            x0 = np.vstack([t ** p for p in range(k)]).ravel()
        if not np.isfinite(f(x0)):
            continue
        used += 1
        r = minimize(f, x0, method="Nelder-Mead", options={"maxiter": 20000, "xatol": 1e-10, "fatol": 1e-12})
        best = min(best, r.fun)
    out.append("k=%d n=%d  V0 positive=%s  minors/P spread=%.1e  E(V0)=%.10f  P(B)/P(J)=%.10f  starts used=%d of 25  best search=%.10f  best>=E0-1e-9: %s"
               % (k, n, d0.min() > 0, ratio_err, e0, formula, used, best, best >= e0 - 1e-9))
print("\n".join(out))
open(sys.argv[1] if len(sys.argv) > 1 else "ogr_check_out.txt", "w").write("\n".join(out) + "\n")
