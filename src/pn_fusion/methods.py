"""Direct vector implementations with explicit pivoted-LU reuse.

No method uses the known root for stopping. Each full solve includes the terminal
residual evaluation; stationary one-macrocycle measurements deliberately do not.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import math
import warnings
import numpy as np
from scipy.linalg import lu_factor, lu_solve, LinAlgWarning
from .models import CostModel, positive_integer
from .problems import Problem

@dataclass
class Resources:
    F: int=0
    J: int=0
    LU: int=0
    solves: int=0
    matvecs: int=0
    matrix_scales: int=0
    vector_scales: int=0

    def work(self, model: CostModel) -> float:
        n=model.n
        return (self.F*n*model.mu0+self.J*n*n*model.mu1+self.LU*model.L
                +self.solves*model.T+(self.matvecs+self.matrix_scales)*n*n
                +self.vector_scales*n)

@dataclass
class Result:
    method: str
    x: np.ndarray
    cycles: int
    residual: float
    error: float
    converged: bool
    resources: Resources
    message: str=""

class Operations:
    def __init__(self, problem: Problem):
        self.problem=problem
        self.stats=Resources()

    def f(self,x):
        self.stats.F+=1
        return self.problem.F(x)

    def j(self,x):
        self.stats.J+=1
        return self.problem.J(x)

    def factor(self,A):
        self.stats.LU+=1
        with warnings.catch_warnings():
            warnings.simplefilter("error",LinAlgWarning)
            return lu_factor(A,overwrite_a=False,check_finite=False)

    def solve(self,factors,b):
        self.stats.solves+=1
        return lu_solve(factors,b,check_finite=False)

    def mv(self,A,b):
        self.stats.matvecs+=1
        return A@b

    def poly(self,factors,K,b,q):
        """Apply H_q b with q solves and q-1 matvecs; never form an inverse."""
        v=self.solve(factors,b); h=v.copy()
        for _ in range(1,q):
            v=v-self.solve(factors,self.mv(K,v))
            h=h+v
        return h


def pn_step(op,x,fx,N,previous=None):
    u=x if previous is None else x-op.solve(previous,fx)
    factors=op.factor(op.j(u))
    z=x-op.solve(factors,fx)
    for _ in range(1,N):
        z=z-op.solve(factors,op.f(z))
    return z,factors


def fused_step(op,x,fx,N,q,previous=None):
    u=x if previous is None else x-op.poly(previous[0],previous[1],fx,q)
    factors=op.factor(op.j(u))
    z=x-op.solve(factors,fx)
    for _ in range(1,N):
        z=z-op.solve(factors,op.f(z))
    y=z; fy=op.f(y)
    v=y-op.solve(factors,fy)
    K=op.j(v)
    w=y-op.poly(factors,K,fy,q)
    for _ in range(1,N):
        w=w-op.poly(factors,K,op.f(w),q)
    return w,(factors,K)


def solve_problem(problem: Problem, method: str, *, N: int=2, q: int|None=None,
                  tol: float=1e-12, maxit: int=20, initialization: str="D") -> Result:
    """Solve using P, S, Phat, Pcompose, M6 or M8.

    A Pcompose iteration is two ordinary P_N cycles. Residual stopping is only
    tested after the complete macrocycle, as in the fused comparison.
    """
    if method not in {"P","S","Phat","Pcompose","M6","M8"}:
        raise ValueError("Unknown method")
    positive_integer(N,"N")
    if method in {"P","Phat","Pcompose"}: positive_integer(N,"N",2)
    q=N+2 if q is None else positive_integer(q,"q",2)
    positive_integer(maxit,"maxit")
    if initialization not in {"I","D"}: raise ValueError("initialization must be I or D")
    if not math.isfinite(tol) or tol<=0: raise ValueError("tol must be finite and positive")
    op=Operations(problem); x=problem.x0.copy(); previous=None
    if initialization=="I" and method in {"P","Pcompose"}:
        previous=op.factor(op.j(x))
    if initialization=="I" and method=="Phat":
        raise ValueError("Immediate fused startup is not used by the archived timing protocol")
    residual=math.inf; fx=op.f(x); cycles=0
    try:
        for cycles in range(1,maxit+1):
            if method=="P":
                x,previous=pn_step(op,x,fx,N,previous)
            elif method=="S":
                x,_=pn_step(op,x,fx,N,None)
            elif method=="Pcompose":
                y,fac=pn_step(op,x,fx,N,previous)
                x,previous=pn_step(op,y,op.f(y),N,fac)
            elif method=="Phat":
                x,previous=fused_step(op,x,fx,N,q,previous)
            elif method=="M6":
                fac=op.factor(op.j(x)); y=x-op.solve(fac,fx); K=op.j(y)
                d=op.solve(fac,op.f(y)); d3=op.solve(fac,op.mv(K,d))
                op.stats.vector_scales+=1; z=y-2*d+d3
                b=op.solve(fac,op.f(z)); b3=op.solve(fac,op.mv(K,b))
                op.stats.vector_scales+=1; x=z-2*b+b3
            else:
                Jx=op.j(x); fac=op.factor(Jx); s=op.solve(fac,fx)
                op.stats.vector_scales+=2
                y=x-0.5*s; z=x-(2/3)*s
                op.stats.matrix_scales+=1
                fac2=op.factor(Jx-3*op.j(z))
                u=y+op.solve(fac2,fx)
                op.stats.vector_scales+=2
                v=u+2*op.solve(fac2,op.f(u))
                x=v+2*op.solve(fac2,op.f(v))
            fx=op.f(x)
            residual=float(np.linalg.norm(fx))
            if not np.isfinite(residual):
                return Result(method,x,cycles,residual,math.inf,False,op.stats,"Nonfinite residual")
            if residual<=tol:
                return Result(method,x,cycles,residual,float(np.linalg.norm(x-problem.alpha)),True,op.stats)
    except (LinAlgWarning, np.linalg.LinAlgError, FloatingPointError, ValueError) as exc:
        return Result(method,x,cycles,residual,float(np.linalg.norm(x-problem.alpha)),False,op.stats,str(exc))
    return Result(method,x,cycles,residual,float(np.linalg.norm(x-problem.alpha)),False,op.stats,"Maximum cycles reached")


def prepare_stationary(problem: Problem, method: str, N: int=2, q: int=4):
    """Build a valid inherited state outside timing; see the provenance limitation.

    The missing historical stationary driver did not archive its inherited
    matrices. This public driver uses the first delayed macrocycle to construct
    them. This is an explicit fresh benchmark, not a reconstruction of old CPU.
    """
    op=Operations(problem); x=problem.x0.copy(); fx=op.f(x)
    if method=="Phat":
        x,previous=fused_step(op,x,fx,N,q)
    elif method=="Pcompose":
        y,fac=pn_step(op,x,fx,N)
        x,previous=pn_step(op,y,op.f(y),N,fac)
    else:
        raise ValueError("stationary mode requires Phat or Pcompose")
    return x,previous


def stationary_macrocycle(problem: Problem, method: str, state, N: int=2, q: int=4):
    op=Operations(problem); x,previous=state; fx=op.f(x)
    if method=="Phat":
        x,previous=fused_step(op,x,fx,N,q,previous)
    else:
        y,fac=pn_step(op,x,fx,N,previous)
        x,previous=pn_step(op,y,op.f(y),N,fac)
    return Result(method,x,1,math.nan,math.nan,True,op.stats,"Stationary macrocycle; no terminal residual")
