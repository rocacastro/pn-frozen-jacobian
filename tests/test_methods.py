import numpy as np
import pytest
from pn_fusion.models import CostModel
from pn_fusion.problems import dense_problem, hammerstein_problem
from pn_fusion.methods import Operations, solve_problem, prepare_stationary, stationary_macrocycle

@pytest.mark.parametrize("field", ['H1','H2','H3','H4','H5'])
def test_analytic_jacobians(field):
    problem=dense_problem(field,delta=0.3)
    n=problem.n
    numerical=np.column_stack([np.imag(problem.F(problem.x0.astype(complex)+1j*1e-100*np.eye(n)[j]))/1e-100 for j in range(n)])
    assert np.linalg.norm(numerical-problem.J(problem.x0))/np.linalg.norm(numerical) < 5e-15

@pytest.mark.parametrize("q", [2,3,4,6])
def test_polynomial_action_and_counts(q):
    problem=hammerstein_problem();op=Operations(problem)
    A=np.array([[3.,.2],[.1,4.]])
    K=A+np.array([[.01,-.02],[.02,-.01]])
    b=np.array([.2,.7]);fac=op.factor(A)
    actual=op.poly(fac,K,b,q)
    E=np.eye(2)-np.linalg.solve(A,K)
    expected=sum((np.linalg.matrix_power(E,j) for j in range(q)))@np.linalg.solve(A,b)
    assert np.allclose(actual,expected,rtol=1e-13,atol=1e-14)
    assert (op.stats.LU,op.stats.solves,op.stats.matvecs)==(1,q,q-1)

@pytest.mark.parametrize("method,N,q", [('P',3,4),('P',5,4),('S',3,4),('Phat',2,4),('Pcompose',2,4),('M6',2,4),('M8',2,4)])
def test_complete_run_counts(method,N,q):
    result=solve_problem(dense_problem('H1',0.6),method,N=N,q=q)
    r=result.cycles;s=result.resources
    assert result.converged and result.residual<=1e-12
    if method=='P': expected=(r*N+1,r,r,r*(N+1)-1,0)
    elif method=='S': expected=(r*N+1,r,r,r*N,0)
    elif method=='Phat': expected=(2*r*N+1,2*r,r,r*(N+1)*(q+1)-q,r*(N+1)*(q-1)-(q-1))
    elif method=='Pcompose': expected=(2*r*N+1,2*r,2*r,2*r*(N+1)-1,0)
    elif method=='M6': expected=(3*r+1,2*r,r,5*r,2*r)
    else: expected=(3*r+1,2*r,2*r,4*r,0)
    assert (s.F,s.J,s.LU,s.solves,s.matvecs)==expected

@pytest.mark.parametrize("method,expected", [('Phat',(4,2,1,15,9)),('Pcompose',(4,2,2,6,0))])
def test_stationary_counts(method,expected):
    p=dense_problem('G1',n=40)
    r=stationary_macrocycle(p,method,prepare_stationary(p,method))
    s=r.resources
    assert (s.F,s.J,s.LU,s.solves,s.matvecs)==expected

def test_first_delayed_cycle_is_shamanskii():
    p=dense_problem('H5',.1)
    a=solve_problem(p,'P',N=10);b=solve_problem(p,'S',N=10)
    assert a.cycles==b.cycles==1
    assert np.array_equal(a.x,b.x)
    assert vars(a.resources)==vars(b.resources)

def test_immediate_initialization_accounts_for_extra_factorization():
    p=dense_problem('H1',.1);r=solve_problem(p,'P',N=3,initialization='I')
    assert r.resources.J==r.resources.LU==r.cycles+1
    assert r.resources.solves==r.cycles*4

@pytest.mark.parametrize("tol,cycles", [(1e-8,1),(1e-10,2),(1e-12,2)])
def test_hammerstein_cycle_switch(tol,cycles):
    r=solve_problem(hammerstein_problem(),'P',N=7,tol=tol)
    assert r.converged and r.cycles==cycles

@pytest.mark.parametrize("kwargs", [{'method':'bad'},{'method':'P','N':1},{'method':'Phat','q':1},{'method':'P','tol':0}])
def test_invalid_solver_arguments(kwargs):
    with pytest.raises(ValueError): solve_problem(hammerstein_problem(),**kwargs)
