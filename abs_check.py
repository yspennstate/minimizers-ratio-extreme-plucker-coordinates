"""Absolute-value version: minimize max|D_I|/min|D_I| over all real k x n matrices (no sign condition), compare with E(V0)."""
import itertools, math, sys
import numpy as np
from scipy.optimize import minimize
rng = np.random.default_rng(5)
out = []
for k, n in [(3, 6), (3, 7)]:
    subs = list(itertools.combinations(range(n), k))
    P = np.array([np.prod([math.sin(math.pi * (b - a) / n) for a, b in itertools.combinations(I, 2)]) for I in subs])
    E0 = P.max() / P.min()
    def f(x):
        X = x.reshape(k, n); d = np.abs(np.array([np.linalg.det(X[:, list(I)]) for I in subs]))
        return np.inf if d.min() < 1e-12 else d.max() / d.min()
    best = np.inf
    for trial in range(60):
        r = minimize(f, rng.standard_normal(k * n), method="Nelder-Mead", options={"maxiter": 40000, "xatol": 1e-10, "fatol": 1e-12})
        r = minimize(f, r.x, method="Nelder-Mead", options={"maxiter": 40000, "xatol": 1e-11, "fatol": 1e-13})
        best = min(best, r.fun)
    out.append("k=%d n=%d E(V0)=%.8f  best |minor| ratio over 60 generic starts=%.8f  below E(V0): %s" % (k, n, E0, best, best < E0 - 1e-6))
    print(out[-1], flush=True)
open(sys.argv[1] if len(sys.argv) > 1 else "abs_check_out.txt", "w").write("\n".join(out) + "\n")
