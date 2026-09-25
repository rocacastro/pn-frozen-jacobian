"""Canonical dense fields and the eight-node Hammerstein system from v0.41."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Callable
import hashlib
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from .models import DIMENSIONS, positive_integer

Array = np.ndarray

@dataclass
class Problem:
    name: str
    n: int
    x0: Array
    alpha: Array
    F: Callable[[Array], Array]
    J: Callable[[Array], Array]
    matrices: tuple[Array, ...] = ()


def canonical_matrices(n: int) -> tuple[Array, Array, Array]:
    positive_integer(n, "n")
    rng = np.random.Generator(np.random.PCG64(20260923+n))
    RA = rng.uniform(-1, 1, size=(n,n))
    RB = rng.uniform(-1, 1, size=(n,n))
    RC = rng.uniform(-1, 1, size=(n,n))
    scale = math.sqrt(n)
    return 4*np.eye(n)+RA/scale, RB/scale, RC/scale


def matrix_hash(matrix: Array) -> str:
    return hashlib.sha256(np.asarray(matrix, dtype='<f8').tobytes(order='C')).hexdigest()


def dense_problem(name: str, delta: float = 0.1, n: int | None = None) -> Problem:
    field = {"G1": "H1", "G5": "H5"}.get(name, name)
    if field not in DIMENSIONS:
        raise ValueError(f"Unknown dense field: {name}")
    if not math.isfinite(delta) or delta < 0:
        raise ValueError("delta must be finite and nonnegative")
    n = DIMENSIONS[field] if n is None else positive_integer(n,"n")
    A,B,C = canonical_matrices(n)
    j = np.arange(1,n+1,dtype=np.float64)
    alpha = j/10
    d = np.sin(np.sqrt(2)*j)+0.5*np.cos(np.sqrt(3)*j)
    x0 = alpha+delta*d/np.linalg.norm(d)
    # Constant matrix scalings are setup work, never timed function work.
    if field == "H1":
        B08,C03,B16,C09 = 0.8*B,0.3*C,1.6*B,0.9*C
        def F(x):
            e=x-alpha; e2=e*e
            return A@e+B08@e2+C03@(e2*e)
        def J(x):
            e=x-alpha
            return A+B16*e[None,:]+C09*(e*e)[None,:]
    elif field == "H2":
        B05,C04=0.5*B,0.4*C
        def F(x):
            e=x-alpha
            return A@e+B05@(np.exp(e)-1-e)+C04@(np.sin(e)-e)
        def J(x):
            e=x-alpha
            return A+B05*(np.exp(e)-1)[None,:]+C04*(np.cos(e)-1)[None,:]
    elif field == "H3":
        B06,C025,C05=0.6*B,0.25*C,0.5*C
        def F(x):
            e=x-alpha
            return A@e+B06@(np.arctan(e)-e)+C025@(e*e)
        def J(x):
            e=x-alpha
            return A+B06*(1/(1+e*e)-1)[None,:]+C05*e[None,:]
    elif field == "H4":
        B05,C03=0.5*B,0.3*C
        def F(x):
            e=x-alpha
            return A@e+B05@np.log1p(e*e)+C03@(np.cos(e)-1+0.5*e*e)
        def J(x):
            e=x-alpha
            return A+B*(e/(1+e*e))[None,:]+C03*(-np.sin(e)+e)[None,:]
    else:
        B04,C03,C06=0.4*B,0.3*C,0.6*C
        def F(x):
            e=x-alpha
            return A@e+B04@(np.exp(e)-1-e)+C03@np.log1p(e*e)
        def J(x):
            e=x-alpha
            return A+B04*(np.exp(e)-1)[None,:]+C06*(e/(1+e*e))[None,:]
    return Problem(name,n,x0,alpha,F,J,(A,B,C))


def hammerstein_problem() -> Problem:
    """Use exactly n=8, x0=0.5*1 and the Green-kernel Nyström discretization."""
    n=8
    nodes,weights=leggauss(n)
    nodes=(nodes+1)/2; weights=weights/2
    s=nodes[:,None]; t=nodes[None,:]
    A=weights[None,:]*np.minimum(s,t)*(1-np.maximum(s,t))
    A3=3*A
    def F(x):
        return 5*x-5-A@(x*x*x)
    def J(x):
        return 5*np.eye(n)-A3*(x*x)[None,:]
    # Independent reference computation outside every timed region.
    alpha=np.ones(n)
    for _ in range(12):
        alpha=alpha-np.linalg.solve(J(alpha),F(alpha))
    return Problem("Hammerstein",n,np.full(n,0.5),alpha,F,J,(A,))
