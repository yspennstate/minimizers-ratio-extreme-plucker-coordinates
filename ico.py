import itertools, numpy as np
p=(1+5**.5)/2
V=np.array([[0,1,p],[0,1,-p],[1,p,0],[1,-p,0],[p,0,1],[-p,0,1]]).T
d=sorted(abs(np.linalg.det(V[:,list(I)])) for I in itertools.combinations(range(6),3))
print("ico 6 axes |minors|:", np.round(d,6), "ratio", d[-1]/d[0], "phi", p)
