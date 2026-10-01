"""Faithful multiprecision implementations for the PN project.

The preceding implementation draft has been retained and completed here:
object arrays, generic partial-pivoting LU reused by triangular solves, the
original x coordinates, delayed start D, and explicit resource accounting.
No inverse, low-rank formula, exact-root predictor or float64 solve is used.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
from pathlib import Path
import hashlib
import math
import numpy as np
from mpmath import mp

FIELDS = {'H1':20, 'H2':50, 'H3':100, 'H4':150, 'H5':200}


def vec(values):
    return np.asarray(list(values), dtype=object)


def norm(x):
    return mp.sqrt(mp.fsum(v*v for v in x))


def text(x, digits=30):
    return mp.nstr(x, digits, strip_zeros=False)


def pack(x):
    s, m, e, b = x._mpf_
    return [int(s), format(int(m), 'x'), int(e), int(b)]


def unpack(t):
    return mp.make_mpf((int(t[0]), int(t[1],16), int(t[2]), int(t[3])))


def pack_vec(x):
    return [pack(v) for v in x]


def unpack_vec(data):
    return vec(unpack(v) for v in data)


def coefficients(n):
    rng = np.random.Generator(np.random.PCG64(20260923+n))
    RA, RB, RC = (rng.uniform(-1,1,(n,n)) for _ in range(3))
    return 4*np.eye(n)+RA/math.sqrt(n), RB/math.sqrt(n), RC/math.sqrt(n)


class DenseProblem:
    """Exact dyadic lift of frozen A/B/C; analytic constants remain decimal.

    x0 is separately frozen at the working precision and lifted exactly for
    validation. alpha_j=j/10 is constructed at each precision. This protocol
    preserves the coefficient realization, not float64 evaluation roundoff.
    """
    def __init__(self, name, n=None, inputs=None):
        if name not in FIELDS:
            raise ValueError(f'Unknown field {name}')
        self.name, self.n = name, FIELDS[name] if n is None else int(n)
        path = None if inputs is None else Path(inputs)/f'coefficients_n{self.n}.npz'
        if path is not None and path.exists():
            with np.load(path, allow_pickle=False) as dat:
                arrays = [dat[k].copy() for k in ('A','B','C')]
        else:
            arrays = coefficients(self.n)
        self.input_sha256 = hashlib.sha256(b''.join(np.asarray(a,dtype='<f8').tobytes(order='C') for a in arrays)).hexdigest()
        self.A, self.B, self.C = [np.asarray([[mp.mpf(float(v)) for v in row] for row in a],dtype=object) for a in arrays]
        self.alpha = vec(mp.mpf(j)/10 for j in range(1,self.n+1))
        fb,fc,jb,jc = {'H1':('0.8','0.3','1.6','0.9'),
                       'H2':('0.5','0.4','0.5','0.4'),
                       'H3':('0.6','0.25','0.6','0.5'),
                       'H4':('0.5','0.3','1','0.3'),
                       'H5':('0.4','0.3','0.4','0.6')}[name]
        # Scalar matrix precomputations are outside all timers.
        self.Bf, self.Cf = self.B*mp.mpf(fb), self.C*mp.mpf(fc)
        self.Bj, self.Cj = self.B*mp.mpf(jb), self.C*mp.mpf(jc)
        self.half = mp.mpf('0.5')

    def initial(self, delta='0.30'):
        d = vec(mp.sin(mp.sqrt(2)*j)+self.half*mp.cos(mp.sqrt(3)*j) for j in range(1,self.n+1))
        return self.alpha+(mp.mpf(str(delta))/norm(d))*d

    def F(self,x):
        e=x-self.alpha; e2=e*e
        if self.name=='H1': v,w=e2,e2*e
        elif self.name=='H2': v,w=vec(mp.expm1(a)-a for a in e),vec(mp.sin(a)-a for a in e)
        elif self.name=='H3': v,w=vec(mp.atan(a)-a for a in e),e2
        elif self.name=='H4': v,w=vec(mp.log1p(a*a) for a in e),vec(mp.cos(a)-1+self.half*a*a for a in e)
        else: v,w=vec(mp.expm1(a)-a for a in e),vec(mp.log1p(a*a) for a in e)
        return self.A@e+self.Bf@v+self.Cf@w

    def J(self,x):
        e=x-self.alpha
        if self.name=='H1': v,w=e,e*e
        elif self.name=='H2': v,w=vec(mp.expm1(a) for a in e),vec(mp.cos(a)-1 for a in e)
        elif self.name=='H3': v,w=vec(1/(1+a*a)-1 for a in e),e
        elif self.name=='H4': v,w=vec(a/(1+a*a) for a in e),vec(a-mp.sin(a) for a in e)
        else: v,w=vec(mp.expm1(a) for a in e),vec(a/(1+a*a) for a in e)
        return self.A+self.Bj*v[None,:]+self.Cj*w[None,:]


class StressProblem:
    name='R2_stress'; n=2; input_sha256='exact-polynomial-R2-v1'
    def __init__(self): self.alpha=vec([mp.mpf(0),mp.mpf(0)])
    def initial(self,delta=None): return vec([mp.mpf('0.02'),mp.mpf('0.04')])
    def F(self,x):
        u,v=x
        return vec([u+u*u+u*v,v+u*u+u*v+v*v])
    def J(self,x):
        u,v=x
        return np.asarray([[1+2*u+v,u],[2*u+v,1+u+2*v]],dtype=object)


@dataclass
class LU:
    a: np.ndarray
    swaps: list


def factor(a):
    a=np.array(a,dtype=object,copy=True); n=len(a)
    if a.shape!=(n,n) or n<1: raise ValueError('A nonempty square matrix is required')
    swaps=[]
    for k in range(n-1):
        p=k+max(range(n-k),key=lambda j:abs(a[k+j,k]))
        if not mp.isfinite(a[p,k]) or not a[p,k]: raise ArithmeticError('Singular or nonfinite pivot')
        if p!=k:
            a[[k,p],:]=a[[p,k],:]; swaps.append((k,p))
        a[k+1:,k]=a[k+1:,k]/a[k,k]
        a[k+1:,k+1:]=a[k+1:,k+1:]-a[k+1:,k,None]*a[None,k,k+1:]
    if not mp.isfinite(a[-1,-1]) or not a[-1,-1]: raise ArithmeticError('Singular final pivot')
    return LU(a,swaps)


def solve(lu,b):
    a,x=lu.a,np.array(b,dtype=object,copy=True); n=len(x)
    for i,j in lu.swaps: x[i],x[j]=x[j],x[i]
    for i in range(1,n): x[i]-=mp.fsum(a[i,j]*x[j] for j in range(i))
    for i in range(n-1,-1,-1): x[i]=(x[i]-mp.fsum(a[i,j]*x[j] for j in range(i+1,n)))/a[i,i]
    return x


@dataclass
class Counts:
    F:int=0; J:int=0; LU:int=0; solves:int=0; matvecs:int=0; scales:int=0
    def dict(self): return asdict(self)


class Operations:
    def __init__(self,p):
        self.p=p; self.c=Counts(); self.half=mp.mpf('0.5'); self.two_thirds=mp.mpf(2)/3
    def F(self,x): self.c.F+=1; return self.p.F(x)
    def J(self,x): self.c.J+=1; return self.p.J(x)
    def factor(self,a): self.c.LU+=1; return factor(a)
    def solve(self,fac,b): self.c.solves+=1; return solve(fac,b)
    def matvec(self,a,b): self.c.matvecs+=1; return a@b
    def poly(self,fac,K,q,b):
        v=self.solve(fac,b); h=v.copy()
        for _ in range(1,q):
            v=v-self.solve(fac,self.matvec(K,v)); h=h+v
        return h


@dataclass(frozen=True)
class Method:
    family:str; N:int=0; q:int=0
    @property
    def key(self):
        if self.family in ('M6','M8'): return self.family
        if self.family=='PH': return f'Phat{self.N}q{self.q}'
        if self.family=='COMP': return f'P{self.N}x2'
        return f'{self.family}{self.N}'


def method(key):
    if key in ('M6','M8'): return Method(key)
    if key.startswith('Phat'):
        n,q=key[4:].split('q'); m=Method('PH',int(n),int(q))
    elif key.startswith('P') and key.endswith('x2'): m=Method('COMP',int(key[1:-2]))
    elif key[0] in ('P','S'): m=Method(key[0],int(key[1:]))
    else: raise ValueError(f'Unknown method {key}')
    if m.N<1 or (m.family in ('P','COMP','PH') and m.N<2) or (m.family=='PH' and m.q<2):
        raise ValueError('Invalid method parameters')
    return m


def pn_step(op,x,Fx,N,memory=None):
    u=x if memory is None else x-op.solve(memory,Fx)
    fac=op.factor(op.J(u)); z=x-op.solve(fac,Fx)
    for _ in range(1,N): z=z-op.solve(fac,op.F(z))
    return z,fac


def step(op,x,Fx,m,memory=None):
    if m.family=='P': return pn_step(op,x,Fx,m.N,memory)
    if m.family=='S': return pn_step(op,x,Fx,m.N,None)[0],None
    if m.family=='COMP':
        y,f=pn_step(op,x,Fx,m.N,memory)
        return pn_step(op,y,op.F(y),m.N,f)
    if m.family=='M6':
        f=op.factor(op.J(x)); y=x-op.solve(f,Fx); K=op.J(y)
        d=op.solve(f,op.F(y)); z=y-2*d+op.solve(f,op.matvec(K,d))
        d=op.solve(f,op.F(z)); op.c.scales+=2*len(x)
        return z-2*d+op.solve(f,op.matvec(K,d)),None
    if m.family=='M8':
        A=op.J(x); f=op.factor(A); s=op.solve(f,Fx)
        y=x-op.half*s; z=x-op.two_thirds*s; K=op.J(z); f2=op.factor(A-3*K)
        u=y+op.solve(f2,Fx); v=u+2*op.solve(f2,op.F(u)); op.c.scales+=len(x)**2+4*len(x)
        return v+2*op.solve(f2,op.F(v)),None
    if m.family=='PH':
        u=x if memory is None else x-op.poly(memory[0],memory[1],m.q,Fx)
        f=op.factor(op.J(u)); z=x-op.solve(f,Fx)
        for _ in range(1,m.N): z=z-op.solve(f,op.F(z))
        y,Fy=z,op.F(z); v=y-op.solve(f,Fy); K=op.J(v)
        w=y-op.poly(f,K,m.q,Fy)
        for _ in range(1,m.N): w=w-op.poly(f,K,m.q,op.F(w))
        return w,(f,K)
    raise ValueError(m)


def expected(m,r,n):
    f=m.family; N=m.N; q=m.q
    if f=='P': out=[r*N+1,r,r,r*(N+1)-1,0,0]
    elif f=='S': out=[r*N+1,r,r,r*N,0,0]
    elif f=='COMP': out=[2*r*N+1,2*r,2*r,2*r*(N+1)-1,0,0]
    elif f=='PH': out=[2*r*N+1,2*r,r,r*(N+1)*(q+1)-q,(r*(N+1)-1)*(q-1),0]
    elif f=='M6': out=[3*r+1,2*r,r,5*r,2*r,2*n*r]
    elif f=='M8': out=[3*r+1,2*r,2*r,4*r,0,r*(n*n+4*n)]
    else: raise ValueError(m)
    return dict(zip(('F','J','LU','solves','matvecs','scales'),out))


def iterate(p,x0,m,tol=None,max_cycles=20,fixed=None,trace=False,retain_state=False):
    """Delayed D; F at the terminal point is cached for the next cycle.

    Only the residual decides stopping. Trace/error/COC are not in the timed
    mode. The initial copy, counters, residual tests and terminal F are timed.
    """
    if max_cycles<1 or (fixed is not None and fixed<1): raise ValueError('At least one cycle required')
    if fixed is None and (tol is None or tol<=0): raise ValueError('Positive tolerance required')
    op=Operations(p); x=np.array(x0,dtype=object,copy=True); memory=None; Fx=op.F(x)
    residuals=[norm(Fx)]; history=[x.copy()] if trace else None
    if not residuals[0]: raise ValueError('Initial point is already a computed root')
    active=0; predictors=0; status='max_cycles'; limit=fixed or max_cycles
    for k in range(limit):
        if m.family in ('P','PH','COMP') and memory is not None: active+=1; predictors+=1
        if m.family=='COMP': predictors+=1
        x,memory=step(op,x,Fx,m,memory); Fx=op.F(x); r=norm(Fx); residuals.append(r)
        if trace: history.append(x.copy())
        if not mp.isfinite(r): status='nonfinite'; break
        if fixed is None and r<=tol: status='converged'; break
        if fixed is not None and not r: status='fixed_computed_zero'; break
    else:
        if fixed is not None: status='fixed_cycles'
    ans={'method':m.key,'status':status,'cycles':k+1,'active_outer_cycles':active,
         'active_predictors':predictors,'counts':op.c.dict(),'x':x,'residual':residuals[-1],
         'residuals':residuals,'history':history}
    if ans['counts']!=expected(m,k+1,p.n): raise AssertionError(f'Resource mismatch: {ans}')
    if retain_state: ans['state']=memory; ans['Fx']=Fx
    return ans


def traces(p,res):
    ee=[norm(x-p.alpha) for x in res['history']]; rows=[]
    for k,x in enumerate(res['history']):
        coc=None
        if k>=2 and all(ee[k-j] for j in (0,1,2)):
            den=mp.log(ee[k-1]/ee[k-2])
            if den: coc=mp.log(ee[k]/ee[k-1])/den
        rows.append({'iteration':k,'error':text(ee[k]),'residual':text(res['residuals'][k]),
                     'coc':None if coc is None else text(coc),
                     'stationary_coc':bool(k>=3 and coc is not None),'x_dyadic':pack_vec(x)})
    return rows
