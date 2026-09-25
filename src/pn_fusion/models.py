"""Algebraic work, reference orders and parameter selection (manuscript Sections 4--6).

q is the number of Neumann terms; the polynomial degree is q - 1.
Reference orders are lower bounds, not assertions about every trajectory.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
from scipy.optimize import brentq

DIMENSIONS = {"H1": 20, "H2": 50, "H3": 100, "H4": 150, "H5": 200}
KAPPA_POINTS = (1.0, 1.2, 1.298, 1.7, 3.23511, 4.45444)


def positive_integer(value: int, name: str, minimum: int = 1) -> int:
    if isinstance(value, bool) or int(value) != value or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def rho(N: int | float) -> float:
    if N < 1 or not math.isfinite(N):
        raise ValueError("N must be finite and >= 1")
    return (N + 2 + math.sqrt((N + 2)**2 - 8)) / 2


def fused_order(N: int, q: int) -> float:
    positive_integer(N, "N", 2)
    positive_integer(q, "q", 2)
    if q <= N:
        trace = (N + 1) * (q + 1)
        return (trace + math.sqrt(trace**2 - 8*q)) / 2
    if q == N + 1:
        trace = (N + 1) * (N + 2)
        return (trace + math.sqrt(trace**2 - 16*(N + 1))) / 2
    return rho(N)**2


def lu_work(n: int, kappa: float = 1.0) -> float:
    positive_integer(n, "n")
    if not math.isfinite(kappa) or kappa < 1:
        raise ValueError("kappa must be finite and >= 1")
    return n*(2*n-1)*(n-1)/6 + kappa*n*(n-1)/2


def solve_work(n: int, kappa: float = 1.0) -> float:
    lu_work(n, kappa)  # Validate the shared domain.
    return n*(n-1) + kappa*n


def field_costs(field: str, n: int, kappa: float = 1.0) -> tuple[float, float]:
    lu_work(n, kappa)
    field = {"G1": "H1", "G5": "H5"}.get(field, field)
    costs = {
        "H1": (3*n+2, 2+1/n), "H2": (3*n+30, 2+30/n),
        "H3": (3*n+15, 2+(1+kappa)/n),
        "H4": (3*n+30, 2+(15+kappa)/n),
        "H5": (3*n+31, 2+(17+kappa)/n),
        "Hammerstein": (n+3, 1+1/n),
    }
    if field not in costs:
        raise ValueError(f"Unknown field: {field}")
    return tuple(float(v) for v in costs[field])


@dataclass(frozen=True)
class CostModel:
    n: int
    mu0: float
    mu1: float
    kappa: float = 1.0

    def __post_init__(self) -> None:
        lu_work(self.n, self.kappa)
        if not all(math.isfinite(v) and v > 0 for v in (self.mu0, self.mu1)):
            raise ValueError("mu0 and mu1 must be finite and positive")

    @classmethod
    def for_field(cls, field: str, n: int | None = None, kappa: float = 1.0):
        if n is None:
            n = 8 if field == "Hammerstein" else DIMENSIONS[field]
        return cls(n, *field_costs(field, n, kappa), kappa)

    @property
    def L(self) -> float:
        return lu_work(self.n, self.kappa)

    @property
    def T(self) -> float:
        return solve_work(self.n, self.kappa)

    def pn(self, N: float) -> float:
        return N*self.n*self.mu0 + self.n**2*self.mu1 + self.L + (N+1)*self.T

    def shamanskii(self, m: float) -> float:
        return m*self.n*self.mu0 + self.n**2*self.mu1 + self.L + m*self.T

    def fused(self, N: int, q: int | None = None) -> float:
        q = N+2 if q is None else q
        return (2*N*self.n*self.mu0 + 2*self.n**2*self.mu1 + self.L
                + (N+1)*(q+1)*self.T + (N+1)*(q-1)*self.n**2)

    def m6(self) -> float:
        return 3*self.n*self.mu0 + 2*self.n**2*self.mu1 + self.L + 5*self.T + 2*self.n**2 + 2*self.n

    def m8(self) -> float:
        return 3*self.n*self.mu0 + 2*self.n**2*self.mu1 + 2*self.L + 4*self.T + self.n**2 + 4*self.n

    def eta_p(self, N: float) -> float:
        return math.log(rho(N)) / self.pn(N)

    def eta_s(self, m: float) -> float:
        return math.log(m+1) / self.shamanskii(m)

    def eta_fused(self, N: int, q: int) -> float:
        return math.log(fused_order(N, q)) / self.fused(N, q)

    def optimal(self, family: str) -> tuple[int, float]:
        """Find the global integer optimum via the proven continuous unimodality."""
        if family == "P":
            lower = 2
            target = (self.n**2*self.mu1+self.L+self.T)/(self.n*self.mu0+self.T)
            fn = lambda x: math.sqrt((x+2)**2-8)*math.log(rho(x))-x-target
            eta = self.eta_p
        elif family == "S":
            lower = 1
            target = (self.n**2*self.mu1+self.L)/(self.n*self.mu0+self.T)
            fn = lambda x: (x+1)*math.log(x+1)-x-target
            eta = self.eta_s
        else:
            raise ValueError("family must be P or S")
        if fn(lower) >= 0:
            return lower, float(lower)
        upper = 2*lower
        while fn(upper) < 0:
            upper *= 2
        x_star = brentq(fn, lower, upper, xtol=1e-13)
        candidates = {max(lower, math.floor(x_star)), math.ceil(x_star)}
        index = max(sorted(candidates), key=eta)
        return index, x_star

    def optimal_q(self, N: int) -> int:
        return max(range(2, N+3), key=lambda q: self.eta_fused(N, q))

    def near_optimal(self, family: str, epsilon: float = 0.005) -> list[int]:
        if not 0 < epsilon < 1:
            raise ValueError("epsilon must lie strictly between zero and one")
        best, _ = self.optimal(family)
        eta = self.eta_p if family == "P" else self.eta_s
        lower = 2 if family == "P" else 1
        threshold = (1-epsilon)*eta(best)
        result = []
        index = lower
        while index <= best or eta(index) >= threshold:
            if eta(index) >= threshold:
                result.append(index)
            index += 1
        return result


def fusion_threshold(N: int, kappa: float = 1.0, delayed: bool = False) -> int:
    """First integer dimension with strictly lower fused algebraic work."""
    positive_integer(N, "N", 2)
    n = 1
    while True:
        difference = (N+1)**2 * (solve_work(n, kappa)+n*n) - lu_work(n, kappa)
        if delayed:
            difference -= (N+1)*(solve_work(n, kappa)+n*n)
        if difference < 0:
            return n
        n += 1


def boundary_p5_m6(n: int, kappa: float = 1.0) -> tuple[float, float, float]:
    """Coefficients of C_P5/log(rho_5) - C_M6/log(6)."""
    lp, lm = math.log(rho(5)), math.log(6)
    return (5*n/lp - 3*n/lm, n*n/lp - 2*n*n/lm,
            (lu_work(n,kappa)+6*solve_work(n,kappa))/lp
            -(lu_work(n,kappa)+5*solve_work(n,kappa)+2*n*n+2*n)/lm)
