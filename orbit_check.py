"""Pointwise test of the per-orbit inequality G_I(V)/G_J(V) >= P(I)/P(J) for every k-subset I,
at random points of Gr>0(k,n) spread out by products of random elementary totally positive bidiagonal factors."""
import itertools, math, sys
import numpy as np
rng = np.random.default_rng(7)

def v0(k, n):
    rows = []
    for r in range(k):
        w = np.exp(1j * math.pi * (k - 1 - 2 * r) / n)
        v = w ** np.arange(n)
        if w.imag > 1e-12: rows += [v.real, v.imag]
        elif abs(w.imag) < 1e-12: rows.append(v.real)
    return np.array(rows[:k])

def tp_factor(n):
    T = np.eye(n)
    for _ in range(3):
        for i in list(range(n - 1)) + list(range(n - 2, -1, -1)):
            for up in (True, False):
                E = np.eye(n); t = math.exp(rng.uniform(-4, 4))
                if up: E[i, i + 1] = t
                else: E[i + 1, i] = t
                T = T @ E
    return T @ np.diag(np.exp(rng.uniform(-3, 3, n)))

out = []
for k, n in [(3, 6), (3, 7), (3, 8), (4, 8), (3, 9), (4, 9), (5, 10)]:
    subs = list(itertools.combinations(range(n), k)); idx = {I: m for m, I in enumerate(subs)}
    orbit = [[idx[tuple(sorted((i + j) % n for i in I))] for j in range(n)] for I in subs]
    P = np.array([np.prod([math.sin(math.pi * (b - a) / n) for a, b in itertools.combinations(I, 2)]) for I in subs])
    M0 = v0(k, n); sgn = np.sign(np.linalg.det(M0[:, :k])); M0[0] *= sgn
    J = idx[tuple(range(k))]
    worst = np.inf; npts = 0; nonpos = 0
    for trial in range(2000):
        X = M0 @ tp_factor(n)
        X = X / np.abs(X).max(axis=1, keepdims=True)
        d = np.array([np.linalg.det(X[:, list(I)]) for I in subs])
        if d.min() <= 0: nonpos += 1; continue
        L = np.log(d); G = np.array([L[o].mean() for o in orbit])
        r = (G - G[J]) - (np.log(P) - np.log(P[J]))
        worst = min(worst, r.min()); npts += 1
    out.append("k=%d n=%d points=%d nonpositive=%d  min over points and I of log[(G_I/G_J)/(P(I)/P(J))] = %.3e  holds=%s"
               % (k, n, npts, nonpos, worst, worst >= -1e-9))
    print(out[-1], flush=True)
open(sys.argv[1] if len(sys.argv) > 1 else "orbit_check_out.txt", "w").write("\n".join(out) + "\n")
