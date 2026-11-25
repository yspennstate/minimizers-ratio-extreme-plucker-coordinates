"""Rank two: the alternating column scaling rho^{(-1)^i/2} keeps E minimal exactly for cos(a)<=rho<=1/cos(a) when n=2 mod 4,
and for no rho != 1 when 4 | n."""
import math, itertools
import numpy as np
def E(n, rho):
    a = math.pi / n
    V = np.array([[math.cos(i*a), math.sin(i*a)] for i in range(n)]).T * np.array([rho**(((-1)**i)/2) for i in range(n)])
    d = [np.linalg.det(V[:, [i, j]]) for i, j in itertools.combinations(range(n), 2)]
    return max(d) / min(d), min(d) > 0
for n in (6, 8, 10, 12, 14, 16):
    a = math.pi / n; e0 = E(n, 1.0)[0]
    lo, hi = math.cos(a), 1 / math.cos(a)
    rows = []
    for rho in (lo*0.999, lo, (lo+1)/2, 1.0, (1+hi)/2, hi, hi*1.001, 1.01, 0.99):
        e, pos = E(n, rho)
        rows.append("%.5f:%s" % (rho, "=" if abs(e-e0) < 1e-9*e0 else ">" if e > e0 else "<"))
    print("n=%d E0=%.8f 1/sin(a)=%.8f  " % (n, e0, 1/math.sin(a)), " ".join(rows))
