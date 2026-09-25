import math
import numpy as np
import pytest
from pn_fusion.models import CostModel, rho, fused_order, lu_work, solve_work, fusion_threshold, boundary_p5_m6

@pytest.mark.parametrize("N", range(2,11))
def test_reference_order_algebra(N):
    p = rho(N)
    assert N+1 < p < N+2
    assert p*p-(N+2)*p+2 == pytest.approx(0, abs=1e-12)
    assert fused_order(N,N+2) == pytest.approx(p*p)
    assert fused_order(N,N+1) < p*p
    assert fused_order(N,N+3) == pytest.approx(p*p)

@pytest.mark.parametrize("n", [1,8,20,55,100,200])
def test_algebraic_work(n):
    assert lu_work(n) == pytest.approx((n**3-n)/3)
    assert solve_work(n) == n*n
    model = CostModel(n,11,2)
    for N in (2,3,6):
        assert model.pn(N)-model.shamanskii(N) == pytest.approx(model.T)
        assert model.fused(N)-2*model.pn(N) == pytest.approx((N+1)**2*(model.T+n*n)-model.L)

@pytest.mark.parametrize("N", [2,3,4,5,6])
def test_integer_threshold(N):
    assert fusion_threshold(N) == 6*(N+1)**2+1
    assert fusion_threshold(N,delayed=True) == 6*N*(N+1)+1

@pytest.mark.parametrize("field,p,s", [('H1',3,3),('H5',9,10)])
def test_independent_optima(field,p,s):
    for kappa in (1,1.2,3.23511,4.5):
        model=CostModel.for_field(field,kappa=kappa)
        assert model.optimal('P')[0] == p
        assert model.optimal('S')[0] == s
        assert model.eta_p(p) >= max(model.eta_p(i) for i in range(2,80))
        assert model.eta_s(s) >= max(model.eta_s(i) for i in range(1,80))

@pytest.mark.parametrize("field,sign", [('H1',1),('H2',1),('H3',-1),('H4',-1),('H5',-1)])
def test_efficiency_boundary(field,sign):
    m=CostModel.for_field(field);a,b,c=boundary_p5_m6(m.n)
    delta=a*m.mu0+b*m.mu1+c
    assert sign*delta>0
    assert delta == pytest.approx(m.pn(5)/math.log(rho(5))-m.m6()/math.log(6))

@pytest.mark.parametrize("args", [(0,1),(8,0),(8,float('nan'))])
def test_invalid_cost_inputs(args):
    with pytest.raises(ValueError): lu_work(*args)

def test_invalid_neumann_terms():
    with pytest.raises(ValueError): fused_order(2,1)
