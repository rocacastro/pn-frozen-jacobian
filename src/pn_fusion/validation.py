"""Integrity, algebraic, and cross-platform deterministic checks.

Timing columns are never used as reproducibility assertions. Small residuals may
round to zero on one platform and remain nonzero on another.
"""
from __future__ import annotations
import csv
import decimal
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from .io import ROOT, read_csv, write_json
from .models import CostModel, DIMENSIONS, rho, fused_order, fusion_threshold, field_costs
from .problems import dense_problem, hammerstein_problem, matrix_hash
from .methods import solve_problem, prepare_stationary, stationary_macrocycle


def require(condition: bool, message: str) -> None:
    if not condition: raise AssertionError(message)


def near(a,b,rtol=2e-11,atol=1e-13):
    return math.isclose(float(a),float(b),rel_tol=rtol,abs_tol=atol)


def check_data_integrity(root: Path=ROOT) -> dict:
    manifest=json.loads((root/'metadata/data_provenance.json').read_text())
    cells=0
    for entry in manifest:
        p=root/entry['path']
        require(hashlib.sha256(p.read_bytes()).hexdigest()==entry['sha256'],f"Modified reference file: {p.name}")
        with p.open(encoding='utf-8',newline='') as f: rows=list(csv.reader(f))
        values=[]
        for i,row in enumerate(rows[1:],1):
            for j,cell in enumerate(row):
                try: decimal.Decimal(cell)
                except decimal.InvalidOperation: continue
                values.append([i,j,cell])
        digest=hashlib.sha256(json.dumps(values,separators=(',',':')).encode()).hexdigest()
        require(digest==entry['numeric_digest'],f"Numeric content changed: {p.name}")
        cells+=len(values)
    return {'datasets':len(manifest),'unaltered_numeric_cells':cells,'passed':True}


def check_models() -> dict:
    checks=0
    for name,N,M in [('p5_m6_indices.csv',5,6),('p7_m8_indices.csv',7,8)]:
        for row in read_csv(name):
            model=CostModel.for_field(row['system'],int(row['n']),float(row['kappa']))
            for value,computed in [(row[f'C_P{N}'],model.pn(N)),(row[f'C_M{M}'],model.m6() if M==6 else model.m8())]:
                require(near(value,computed),f"Cost mismatch: {name}, {row}")
                checks+=1
    for row in read_csv('phat2_p10_stationary_model.csv'):
        model=CostModel.for_field(row['system'],int(row['n']),float(row['kappa']))
        for key,value in [('cost_Phat2',model.fused(2)),('cost_P10',model.pn(10))]:
            require(near(row[key],value),f"Model mismatch: {key}");checks+=1
    for row in read_csv('hammerstein_stationary_model.csv'):
        model=CostModel.for_field('Hammerstein',8,float(row['kappa']))
        for key,value in [('cost_P5',model.pn(5)),('cost_P7',model.pn(7)),('cost_M6',model.m6()),('cost_M8',model.m8())]:
            require(near(row[key],value),f"Hammerstein cost mismatch: {key}");checks+=1
    for row in read_csv('q_sweep.csv'):
        model=CostModel.for_field(row['system'],int(row['n']),float(row['kappa']))
        N,q=int(row['N']),int(row['q'])
        for key,value in [('p_reference',fused_order(N,q)),('cost',model.fused(N,q)),('eta',model.eta_fused(N,q))]:
            require(near(row[key],value),f"q-sweep mismatch: {key}");checks+=1
    for row in read_csv('kappa_pn_shamanskii.csv'):
        model=CostModel.for_field(row['system'],kappa=float(row['kappa']))
        require(model.optimal('P')[0]==int(row['N_P_star']),'P optimum mismatch')
        require(model.optimal('S')[0]==int(row['m_S_star']),'S optimum mismatch');checks+=2
    for row in read_csv('kappa_fusion_thresholds.csv'):
        if row['family'].startswith('Phat'):
            require(fusion_threshold(int(row['index']),float(row['kappa']))==int(row['n_min']),'Fusion threshold mismatch');checks+=1
    return {'scalar_checks':checks,'passed':True}


def check_fingerprints() -> dict:
    matches=0
    for row in read_csv('matrix_fingerprints.csv'):
        p=dense_problem(row['system'],n=int(row['n']))
        for name,A in zip(('A','B','C'),p.matrices):
            require(matrix_hash(A)[:16]==row['sha16_'+name],f"Matrix fingerprint mismatch: {row['system']} {name}")
            matches+=1
    return {'matrix_fingerprints':matches,'passed':True}


def check_trajectories() -> dict:
    checked=0; differences=[]; cache={}
    def get_problem(field,n,delta):
        key=(field,n,round(delta,12))
        if key not in cache:
            cache[key]=hammerstein_problem() if field=='Hammerstein' else dense_problem(field,delta,n)
        return cache[key]
    datasets=('pn_shamanskii_results.csv','p5_m6_results.csv','p7_m8_results.csv',
              'fusion_end_to_end_results.csv','hammerstein_results.csv','phat2_p10_results.csv')
    with threadpool_limits(limits=1,user_api='blas'):
        for name in datasets:
            for row in read_csv(name):
                field=row.get('system',row.get('family','Hammerstein'))
                n=int(row.get('n','8')); delta=float(row.get('delta','0.1'))
                if name=='pn_shamanskii_results.csv':
                    field=row['system']; method=row['family']; N=int(row['index'])
                else:
                    raw=row['method']; N=2
                    if raw in ('P5','P7','P10'): method='P';N=int(raw[1:])
                    elif raw in ('hatP2','Phat2'): method='Phat'
                    elif raw in ('P2circP2','P2oP2','P2_P2','P2composeP2'): method='Pcompose'
                    else: method=raw
                tol=float(row.get('tol','1e-12'))
                problem=get_problem(field,n,delta)
                result=solve_problem(problem,method,N=N,tol=tol)
                expected_cycles=int(row.get('cycles',row.get('macrocycles','0')))
                if result.cycles!=expected_cycles or result.converged!=bool(int(row['converged'])):
                    differences.append(dict(dataset=name,field=field,n=n,delta=delta,tol=tol,method=method,N=N,
                                            archived_cycles=expected_cycles,current_cycles=result.cycles))
                require(result.converged and result.residual<=tol,f"Residual criterion failed: {name} {row}")
                for key,actual in [('F_total_including_terminal',result.resources.F),('J',result.resources.J),
                                   ('LU',result.resources.LU),('solves',result.resources.solves),('matvecs',result.resources.matvecs)]:
                    if key in row:
                        require(int(row[key])==actual,f"Resource mismatch {key}: {name} {row} vs {actual}")
                if 'matvec_entries' in row: require(int(row['matvec_entries'])==result.resources.matvecs*n*n,'Matvec entries mismatch')
                if 'matrix_scale_entries' in row: require(int(row['matrix_scale_entries'])==result.resources.matrix_scales*n*n,'Matrix scaling mismatch')
                if 'vector_scale_entries' in row: require(int(row['vector_scale_entries'])==result.resources.vector_scales*n,'Vector scaling mismatch')
                checked+=1
    require(not differences,f"Cycle differences detected: {differences}")
    return {'archived_configurations_checked':checked,'cycle_mismatches':differences,'passed':True,
            'scope':'Reimplemented solvers; residual criterion and resource counts, not bitwise residuals or historical timings'}


def check_stationary_resources() -> dict:
    checked=0
    for row in read_csv('fusion_stationary_results.csv'):
        method='Phat' if row['method']=='hatP2' else 'Pcompose'
        problem=dense_problem(row['family'],0.1,int(row['n']))
        result=stationary_macrocycle(problem,method,prepare_stationary(problem,method))
        for key in ('F','J','LU','solves','matvecs'):
            require(int(row[key])==getattr(result.resources,key),'Stationary resource mismatch')
        checked+=1
    return {'configurations':checked,'passed':True,'historical_state_reconstructed':False}


def run_validation(output: Path, trajectories: bool=True) -> dict:
    report={'reference_integrity':check_data_integrity(),'model':check_models(),'matrices':check_fingerprints()}
    if trajectories:
        report['trajectories']=check_trajectories()
        report['stationary_resources']=check_stationary_resources()
    write_json(output,report)
    return report


def audit_archived_work(output: Path) -> dict:
    """Compare archived work with the current formulas without rewriting any data."""
    from .methods import Resources
    from .io import write_csv
    rows=[];checked=0
    for name in ('pn_shamanskii_results.csv','p5_m6_results.csv','p7_m8_results.csv',
                 'fusion_stationary_results.csv','fusion_end_to_end_results.csv',
                 'hammerstein_results.csv','phat2_p10_results.csv'):
        for row_index,row in enumerate(read_csv(name)):
            n=int(row.get('n','8'));field=row.get('system',row.get('family','Hammerstein'))
            s=Resources(F=int(row.get('F_total_including_terminal',row.get('F','0'))),J=int(row['J']),LU=int(row['LU']),solves=int(row['solves']),
                        matvecs=int(row.get('matvecs','0'))+int(row.get('matvec_entries','0'))//(n*n),
                        matrix_scales=int(row.get('matrix_scales','0'))+int(row.get('matrix_scale_entries','0'))//(n*n),
                        vector_scales=int(row.get('vector_scales','0'))+int(row.get('vector_scale_entries','0'))//n)
            for key in row:
                if key.startswith('work_total_k'): kappa=float(key[len('work_total_k'):].replace('p','.'))
                elif key.startswith('work_k1p') or key.startswith('work_k3p'): kappa=float(key[len('work_k'):].replace('p','.'))
                elif key=='work_kappa1': kappa=1.0
                else: continue
                current=s.work(CostModel.for_field(field,n,kappa));archived=float(row[key]);difference=archived-current
                checked+=1
                if abs(difference)>1e-7:
                    rows.append(dict(dataset=name,row_index=row_index,system=field,n=n,kappa=kappa,method=row.get('method',row.get('family','')),
                                     archived_work=archived,canonical_work=current,archived_minus_canonical=difference,
                                     relative_discrepancy_pct=100*difference/current,
                                     explanation='Archived G5 mu1 fixed at 2+18.2/n; v0.41 declares 2+(17+kappa)/n' if field=='G5' else 'Needs investigation'))
    if rows: write_csv(Path(output)/'archived_work_discrepancies.csv',rows)
    summary={'scalar_work_checks':checked,'discrepant_values':len(rows),
             'affected_datasets':sorted({r['dataset'] for r in rows}),
             'max_absolute_relative_discrepancy_pct':max([abs(r['relative_discrepancy_pct']) for r in rows],default=0),
             'archived_data_modified':False,'status':'Known discrepancy documented' if rows else 'All checks passed'}
    write_json(Path(output)/'archived_work_audit.json',summary)
    return summary


def check_high_precision_output(directory: Path) -> dict:
    comparisons=[]
    for name,column,tolerance in [('high_precision_orders.csv','final_COC',2e-11),('q_order_verification.csv','COC_observed',1e-8)]:
        with (Path(directory)/name).open(newline='') as f: computed=list(csv.DictReader(f))
        archived=read_csv(name)
        require(len(computed)==len(archived),'High-precision row count mismatch')
        for a,b in zip(archived,computed):
            difference=abs(float(a[column])-float(b[column]))
            require(difference<tolerance,'High-precision COC differs beyond archived rounding')
            comparisons.append(dict(dataset=name,N=int(a['N']),q=int(a.get('q',0)),absolute_difference=difference,tolerance=tolerance))
    report={'checked_rows':len(comparisons),'passed':True,'comparisons':comparisons}
    write_json(Path(directory)/'verification.json',report)
    return report
