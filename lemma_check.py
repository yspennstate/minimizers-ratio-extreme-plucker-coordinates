import itertools, random, numpy as np
random.seed(3)
np.random.seed(3)
def tp(k,n):
    t=np.sort(np.random.rand(n))+np.arange(n)
    return np.vstack([t**p for p in range(k)])
def D(X,I): return np.linalg.det(X[:,sorted(I)])
bad=0
for k,n in [(3,7),(4,8),(4,9),(5,10)]:
    X=tp(k,n)
    for _ in range(300):
        a,b,c,d=sorted(random.sample(range(n),4)); rest=[i for i in range(n) if i not in (a,b,c,d)]
        S=random.sample(rest,k-2)
        l=D(X,S+[a,c])*D(X,S+[b,d]); r=D(X,S+[a,b])*D(X,S+[c,d])+D(X,S+[a,d])*D(X,S+[b,c])
        if abs(l-r)>1e-8*abs(l): bad+=1
print("three-term failures:",bad)
from math import gcd
ok=True
for n in range(3,25):
    for k in range(1,n):
        B={(j*n)//k for j in range(k)}
        for r in range(n):
            T={(r+j*n)//k for j in range(k)}
            if not any(T=={(x+e)%n for x in B} for e in range(n)): ok=False; print("not translate",k,n,r)
        # row multiset union check: each element k times
print("T_r translates of B for all 1<=k<n<=24:",ok)
