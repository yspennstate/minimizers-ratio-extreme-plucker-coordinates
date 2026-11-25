"""(a) column scalings of V0: E stays at E(V0) exactly on the subspace x d-periodic with residue-class sums mod g zero,
and rises off it; (b) sorting bound G_I <= G_B at random positive points."""
import itertools, math, sys
import numpy as np
rng = np.random.default_rng(11)
def v0(k, n):
    rows = []
    for r in range(k):
        w = np.exp(1j * math.pi * (k - 1 - 2 * r) / n); v = w ** np.arange(n)
        if w.imag > 1e-12: rows += [v.real, v.imag]
        elif abs(w.imag) < 1e-12: rows.append(v.real)
    M = np.array(rows[:k])
    if np.linalg.det(M[:, :k]) < 0: M[0] *= -1
    return M
def tp_factor(n):
    T = np.eye(n)
    for _ in range(3):
        for i in list(range(n - 1)) + list(range(n - 2, -1, -1)):
            for up in (True, False):
                E = np.eye(n); t = math.exp(rng.uniform(-4, 4))
                if up: E[i, i + 1] = t
                else: E[i + 1, i] = t
                T = T @ E
    return T
out = []
for k, n in [(2, 6), (2, 8), (2, 10), (3, 6), (3, 9), (3, 12), (4, 8), (4, 12), (6, 12)]:
    d = math.gcd(k, n); g = math.gcd(d, n // d)
    subs = list(itertools.combinations(range(n), k)); idx = {I: m for m, I in enumerate(subs)}
    M0 = v0(k, n); d0 = np.array([np.linalg.det(M0[:, list(I)]) for I in subs]); E0 = d0.max() / d0.min()
    def Ex(x):
        dd = d0 * np.exp(np.array([x[list(I)].sum() for I in subs])); return dd.max() / dd.min()
    # basis of the predicted subspace: d-periodic, sum over each residue class mod g of one period = 0
    basis = []
    for j in range(d):
        for j2 in range(j + 1, d):
            if j % g == j2 % g:
                y = np.zeros(d); y[j] = 1; y[j2] = -1; basis.append(y)
    B = np.array(basis).reshape(-1, d)
    rank = np.linalg.matrix_rank(B) if len(basis) else 0
    on = []; off = []; nonper = []
    for _ in range(40):
        if rank:
            y = rng.standard_normal(B.shape[0]) @ B; y *= 1e-3 / np.abs(y).max()
            on.append(Ex(np.tile(y, n // d)) / E0 - 1)
        z = rng.standard_normal(d); z -= z.mean(); 
        if g > 1 or True:
            # period-d with nonzero class sums mod g (only possible when g > 1)
            if g > 1:
                cls = np.array([z[[j for j in range(d) if j % g == t]].sum() for t in range(g)])
                if np.abs(cls).max() > 1e-6: off.append(Ex(np.tile(1e-3 * z / np.abs(z).max(), n // d)) / E0 - 1)
        w = rng.standard_normal(n); w -= w.mean(); nonper.append(Ex(1e-3 * w / np.abs(w).max()) / E0 - 1)
    # sorting bound
    Bset = tuple(sorted(math.floor(j * n / k) for j in range(k)))
    orbit = lambda I: [idx[tuple(sorted((i + t) % n for i in I))] for t in range(n)]
    ob = orbit(Bset); orbs = [orbit(I) for I in subs]; worst = np.inf; npts = 0
    for _ in range(150):
        X = M0 @ tp_factor(n); dd = np.array([np.linalg.det(X[:, list(I)]) for I in subs])
        if dd.min() <= 0: continue
        npts += 1
        L = np.log(dd); GB = L[ob].mean(); worst = min(worst, min(GB - L[o].mean() for o in orbs))
    out.append("k=%d n=%d d=%d g=%d predicted dim=%d basis rank=%d | on-subspace samples=%d max |rel dE|=%s | period-d off-subspace samples=%d min rel dE=%s | generic samples=%d min rel dE=%.1e | sorting points=%d min log(G_B/G_I)=%.2e"
               % (k, n, d, g, d - g, rank, len(on), ("%.1e" % max(abs(v) for v in on)) if on else "n/a",
                  len(off), ("%.1e" % min(off)) if off else "n/a", len(nonper), min(nonper), npts, worst))
    print(out[-1], flush=True)
open(sys.argv[1] if len(sys.argv) > 1 else "col_check_out.txt", "w").write("\n".join(out) + "\n")
