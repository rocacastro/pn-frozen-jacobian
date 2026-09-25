"""Regenerate parameter studies, interval checks, and structural diagnostics."""
from __future__ import annotations
from pathlib import Path
import math
import numpy as np
from scipy.optimize import brentq
from threadpoolctl import threadpool_limits
from .models import CostModel,DIMENSIONS,KAPPA_POINTS,rho,fused_order,fusion_threshold
from .problems import dense_problem,matrix_hash
from .io import write_csv,write_json


def parameter_studies(output: Path) -> dict:
    output=Path(output)
    orders=[];optima=[];bands=[];curvatures=[];qs=[];thresholds=[];comparators=[]
    for kappa in KAPPA_POINTS:
        for field,n in DIMENSIONS.items():
            model=CostModel.for_field(field,n,kappa)
            comparators.append(dict(system=field,n=n,kappa=kappa,
                eta_P5=model.eta_p(5),eta_M6=math.log(6)/model.m6(),
                eta_P7=model.eta_p(7),eta_M8=math.log(8)/model.m8()))
            for N in range(2,7):
                thresholds.append(dict(kappa=kappa,N=N,n_min_stationary=fusion_threshold(N,kappa),n_min_delayed=fusion_threshold(N,kappa,True))) if field=='H1' else None
                for q in range(2,N+3):
                    qs.append(dict(kappa=kappa,system=field,n=n,N=N,q=q,p_reference=fused_order(N,q),
                                   cost=model.fused(N,q),eta=model.eta_fused(N,q),is_optimal=int(q==model.optimal_q(N))))
            if field not in ('H1','H5'): continue
            NP,xP=model.optimal('P');mS,xS=model.optimal('S')
            optima.append(dict(system=field,kappa=kappa,N_P_star=NP,m_S_star=mS,etaP_star=model.eta_p(NP),
                               etaS_star=model.eta_s(mS),relative_advantage_P_pct=100*(model.eta_p(NP)/model.eta_s(mS)-1)))
            for eps in (0.001,0.0025,0.005,0.01,0.02,0.05):
                bands.append(dict(system=field,kappa=kappa,epsilon=eps,P_band=','.join(map(str,model.near_optimal('P',eps))),
                                  S_band=','.join(map(str,model.near_optimal('S',eps)))))
            D=math.sqrt((xP+2)**2-8)
            chiP=(xP+2)/(D**3*math.log(rho(xP)))
            chiS=1/((xS+1)**2*math.log(xS+1))
            level=.995*model.eta_p(xP)
            left=2.0 if model.eta_p(2)>=level else brentq(lambda x:model.eta_p(x)-level,2,xP)
            hi=xP*2
            while model.eta_p(hi)>level: hi*=2
            right=brentq(lambda x:model.eta_p(x)-level,xP,hi)
            curvatures.append(dict(system=field,kappa=kappa,N_star=NP,x_P_star=xP,m_star=mS,x_S_star=xS,
                                   normalized_curvature_P=chiP,normalized_curvature_S=chiS,
                                   continuous_band_left=left,continuous_band_right=right,continuous_band_width=right-left))
    for name,rows in [('family_optima',optima),('near_optimal_bands',bands),('curvature',curvatures),
                      ('q_sweep',qs),('fusion_thresholds',thresholds),('external_indices',comparators)]:
        write_csv(output/(name+'.csv'),rows)
    report=interval_certificates()
    write_json(output/'interval_checks.json',report)
    return {'optima':len(optima),'q_cases':len(qs),'interval_checks':len(report)}


def interval_certificates() -> list[dict]:
    """Endpoint checks for affine cost differences; not a grid-only argument.

    The global integer optimum follows from unimodality: for a fixed candidate,
    checking its two neighbors throughout the interval is sufficient.
    """
    results=[]
    for field,n in DIMENSIONS.items():
        models=[CostModel.for_field(field,n,k) for k in (1.0,4.5)]
        for N,other,p in [(5,'m6',6),(7,'m8',8)]:
            values=[m.pn(N)/math.log(rho(N))-getattr(m,other)()/math.log(p) for m in models]
            results.append(dict(claim=f'{field}: P{N}/{other.upper()}',endpoint_values=values,constant_sign=values[0]*values[1]>0))
        for N in range(2,7):
            best=models[0].optimal_q(N)
            values=[m.fused(N,best)/math.log(fused_order(N,best))-m.fused(N,q)/math.log(fused_order(N,q))
                    for m in models for q in range(2,N+3) if q!=best]
            results.append(dict(claim=f'{field}: optimal q for N={N}',q_star=best,max_normalized_cost_difference=max(values),passed=max(values)<0))
        if field in ('H1','H5'):
            for family in ('P','S'):
                best=models[0].optimal(family)[0];lower=2 if family=='P' else 1
                values=[]
                for m in models:
                    cost=m.pn if family=='P' else m.shamanskii
                    order=rho if family=='P' else lambda j:j+1
                    for candidate in (best-1,best+1):
                        if candidate>=lower: values.append(cost(best)/math.log(order(best))-cost(candidate)/math.log(order(candidate)))
                results.append(dict(claim=f'{field}: global {family} optimum by neighbor comparison and unimodality',
                                    index=best,max_normalized_cost_difference=max(values),passed=max(values)<0))
    return results


def structural_audit(output: Path) -> dict:
    rows=[];fingerprints=[]
    with threadpool_limits(limits=1,user_api='blas'):
        for field,n in DIMENSIONS.items():
            for delta in (.1,.3,.6):
                p=dense_problem(field,delta,n);J=p.J(p.x0);off=J-np.diag(np.diag(J));nonlinear=J-p.matrices[0]
                def residual_rank3(A):
                    s=np.linalg.svd(A,compute_uv=False)
                    return float(np.linalg.norm(s[3:])/np.linalg.norm(s))
                h=1e-100
                derivative=np.empty_like(J)
                for j in range(n):
                    z=p.x0.astype(complex);z[j]+=1j*h;derivative[:,j]=p.F(z).imag/h
                rows.append(dict(system=field,n=n,delta=delta,rank_J=int(np.linalg.matrix_rank(J)),
                    rank_offdiagonal=int(np.linalg.matrix_rank(off)),rank_nonlinear=int(np.linalg.matrix_rank(nonlinear)),
                    density=float(np.count_nonzero(J)/J.size),condition_2=float(np.linalg.cond(J)),
                    offdiagonal_rank3_residual=residual_rank3(off),nonlinear_rank3_residual=residual_rank3(nonlinear),
                    complex_step_relative_error=float(np.linalg.norm(J-derivative)/np.linalg.norm(J))))
            fingerprints.append(dict(system=field,n=n,seed=20260923+n,**{'sha256_'+key:matrix_hash(A) for key,A in zip(('A','B','C'),p.matrices)}))
    write_csv(Path(output)/'structural_audit.csv',rows);write_csv(Path(output)/'matrix_fingerprints.csv',fingerprints)
    return {'configurations':len(rows),'max_complex_step_relative_error':max(r['complex_step_relative_error'] for r in rows)}
