"""Independent matrix-form comparisons and deterministic correctness tests."""
from __future__ import annotations
import numpy as np
from mpmath import mp
from .core import *
from .costs import optimum, eta
from .experiments import js, env


def reference(p,x0,m,cycles=2):
    """Small-matrix oracle: explicit inverses and powers ONLY in untimed tests.

    This realization does not call the production step, LU or solve routines.
    """
    def F(x): return mp.matrix(list(p.F(vec(x))))
    def J(x): return mp.matrix(p.J(vec(x)).tolist())
    def pn(x,N,Aold):
        u=x if Aold is None else x-Aold**-1*F(x)
        A=J(u); z=mp.matrix(x)
        for _ in range(N): z=z-A**-1*F(z)
        return z,A
    x=mp.matrix(list(x0)); old=None; xs=[mp.matrix(x)]
    for _ in range(cycles):
        if m.family=='P': x,old=pn(x,m.N,old)
        elif m.family=='S': x,_=pn(x,m.N,None)
        elif m.family=='COMP':
            y,A=pn(x,m.N,old); x,old=pn(y,m.N,A)
        elif m.family=='M6':
            A=J(x); y=x-A**-1*F(x); H=(2*mp.eye(p.n)-A**-1*J(y))*A**-1
            z=y-H*F(y); x=z-H*F(z)
        elif m.family=='M8':
            A=J(x); y=x-mp.mpf('0.5')*A**-1*F(x); z=(4*y-x)/3
            B=(A-3*J(z))**-1; u=y+B*F(x); v=u+2*B*F(u); x=v+2*B*F(v)
        elif m.family=='PH':
            u=x if old is None else x-old*F(x); A=J(u); z=mp.matrix(x)
            for _ in range(m.N): z=z-A**-1*F(z)
            y=z; v=y-A**-1*F(y); K=J(v); E=mp.eye(p.n)-A**-1*K
            H=sum((E**j for j in range(m.q)),mp.zeros(p.n))*A**-1
            w=mp.matrix(y)
            for _ in range(m.N): w=w-H*F(w)
            x,old=w,H
        xs.append(mp.matrix(x))
    return xs


def self_test(out):
    passed=[]
    def check(cond,label):
        if not cond: raise AssertionError(label)
        passed.append(label)
    with mp.workdps(160):
        rng=np.random.default_rng(917)
        for n in (2,3,5,7):
            A=np.asarray([[mp.mpf(int(v))/7 for v in row] for row in rng.integers(-8,9,(n,n))],dtype=object)
            A+=mp.mpf(12)*np.eye(n,dtype=object); A[[0,-1]]=A[[-1,0]]
            lu=factor(A)
            for trial in range(3):
                b=vec(mp.mpf(int(v))/11 for v in rng.integers(-9,10,n))
                x=solve(lu,b); ref=mp.lu_solve(mp.matrix(A.tolist()),mp.matrix(list(b)))
                check(norm(A@x-b)<mp.mpf('1e-150'),f'LU residual n={n},rhs={trial}')
                check(norm(x-vec(ref))<mp.mpf('1e-150'),f'LU oracle n={n},rhs={trial}')
        try: factor(np.asarray([[mp.mpf(0),mp.mpf(0)],[mp.mpf(0),mp.mpf(0)]],dtype=object))
        except ArithmeticError: check(True,'Singular matrix rejected')
        else: check(False,'Singular matrix rejected')
        A=mp.matrix([[4,1],[2,5]]); K=mp.matrix([[mp.mpf('4.1'),mp.mpf('1.2')],[mp.mpf('1.8'),mp.mpf('5.1')]])
        check(mp.norm(A*K-K*A)>mp.mpf('.01'),'Noncommuting polynomial test matrices')
        p=StressProblem()
        for q in range(2,7):
            op=Operations(p); f=op.factor(np.asarray(A.tolist(),dtype=object)); b=vec([mp.mpf(1),mp.mpf(2)])
            h=op.poly(f,np.asarray(K.tolist(),dtype=object),q,b); E=mp.eye(2)-A**-1*K
            H=sum((E**j for j in range(q)),mp.zeros(2))*A**-1
            check(norm(h-vec(H*mp.matrix(list(b))))<mp.mpf('1e-150'),f'Polynomial action q={q}')
            check(mp.norm(H*K-(mp.eye(2)-E**q))<mp.mpf('1e-150'),f'Polynomial defect q={q}')
            check(op.c.solves==q and op.c.matvecs==q-1,f'Polynomial resources q={q}')
        for name in FIELDS:
            p=DenseProblem(name,3); x=p.initial('0.30'); J=p.J(x)
            for j in range(3):
                def column(t):
                    z=x.copy(); z[j]=t
                    return p.F(z)
                for i in range(3):
                    d=mp.diff(lambda t:column(t)[i],x[j])
                    check(abs(d-J[i,j])<mp.mpf('1e-145'),f'Analytic Jacobian {name}:{i},{j}')
            check(norm(p.F(p.alpha))==0,f'Known root {name}')
            for key in ('P2','P3','S2','S3','P5','M6','P7','M8','Phat2q2','Phat2q3','Phat2q4','P2x2','P10'):
                m=method(key); actual=iterate(p,x,m,fixed=2,trace=True); rr=reference(p,x,m,2)
                for k,z in enumerate(actual['history']):
                    check(norm(z-vec(rr[k]))<mp.mpf('1e-140'),f'Matrix-form oracle {name}:{key}:k={k}')
                check(actual['counts']==expected(m,actual['cycles'],p.n),f'Full resources {name}:{key}')
            a=iterate(p,x,method('P3'),fixed=1); b=iterate(p,x,method('S3'),fixed=1)
            check(pack_vec(a['x'])==pack_vec(b['x']),f'Identical delayed first block {name}')
        x=vec([mp.mpf('0.02'),mp.mpf('1e-1024'),mp.mpf('-9.3')]); packed=pack_vec(x)
        with mp.workdps(20): restored=unpack_vec(packed)
        check(pack_vec(restored)==packed,'Exact dyadic serialization across precisions')
        for name,n in FIELDS.items():
            for fam,first in [('P',2),('S',1)]:
                i=optimum(fam,name,n); brute=max(range(first,201),key=lambda j:eta(f'{fam}{j}',name,n))
                check(i==brute,f'Independent optimum {name}:{fam}')
        for q,s,v in [(2,9,3),(3,12,6),(4,15,9)]:
            m=method(f'Phat2q{q}'); a,b=expected(m,1,7),expected(m,2,7)
            check(b['solves']-a['solves']==s and b['matvecs']-a['matvecs']==v,f'Fused stationary counts q={q}')
    report={'passed':True,'checks':len(passed),'details':passed,'environment':env(),
            'scope':'Untimed numerical implementation tests, not a proof or completed timing campaign'}
    js(out,report); print(f'{len(passed)} deterministic checks passed.',flush=True)
    return report
