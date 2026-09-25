"""English table exports from immutable data and the declared algebraic model.

Table numbers and source labels refer to manuscript v0.41. Timing exports use the
archived measurements, never the computer running this script.
"""
from __future__ import annotations
from pathlib import Path
import csv
import json
import math
import statistics
from .io import ROOT,read_csv,write_csv,write_json
from .models import CostModel,DIMENSIONS,KAPPA_POINTS,rho,fused_order,fusion_threshold

TABLES=[
('tab:kappa-threshold-sensitivity','Stationary fusion threshold sensitivity'),
('tab:q-transition','Reference orders and the q-transition trade-off'),
('tab:qstar-fields','Optimal number of Neumann terms'),
('tab:precision-regimes','Experimental protocols'),
('tab:H-fields','Canonical dense fields and evaluation costs'),
('tab:elementary-products','Elementary equivalent-product weights'),
('tab:H-structural-audit','Structural audit of the dense fields'),
('tab:PN-initialization-results','Immediate and delayed initialization'),
('tab:hammerstein-stationary','Hammerstein stationary model'),
('tab:hammerstein-practical','Hammerstein end-to-end comparison'),
('tab:optimal-P-S-kappa','Independently optimized families'),
('tab:obj12-flatness','Flatness of the stationary optimum'),
('tab:asymptotic-members-finite','Stationarily optimal members at finite tolerance'),
('tab:practical-minima-P-S','Observed practical minima of the family sweep'),
('tab:external-indices','Stationary indices of the external comparators'),
('tab:external-cpu','External comparisons over all configurations'),
('tab:fusion-resource-accounting','Stationary and transient fusion resources'),
('tab:fusion-stationary','Stationary fusion comparison'),
('tab:fusion-practical','End-to-end fusion comparison'),
('tab:Phat2-P10','Secondary comparison: Phat2 versus P10'),
('tab:order-verification','Independent high-precision order verification'),
('tab:q-coc','Order verification across the q transition')]


def numeric(row,key): return float(row[key])


def table_rows() -> list[tuple[list[dict],list[str]]]:
    result=[]
    rows=[]
    for N in range(2,7):
        row={'N':N}
        for k in KAPPA_POINTS: row[f'kappa_{k:g}']=fusion_threshold(N,k)
        rows.append(row)
    result.append((rows,['models.fusion_threshold']))
    rows=[]
    for N in range(2,7):
        sigma=fused_order(N,N+1);full=rho(N)**2
        rows.append(dict(N=N,rho=rho(N),sigma=sigma,full_order=full,relative_order_gain_pct=100*(full/sigma-1),Xi=math.log(sigma)/math.log(full/sigma)))
    result.append((rows,['models.rho','models.fused_order']))
    result.append(([dict(system=field,**{f'N_{N}':CostModel.for_field(field).optimal_q(N) for N in range(2,7)}) for field in DIMENSIONS],['models.CostModel.optimal_q']))
    rows=[dict(block='Initialization I/D',tolerance='1e-12',precision='float64',warmup=3,repetitions=51,order='alternating'),
          dict(block='PN/S family sweep',tolerance='1e-12',precision='float64',warmup=2,repetitions=31,order='pseudorandom')]
    for block in ('P5/M6','P7/M8','Hammerstein','Fusion','Phat2/P10'):
        rows.append(dict(block=block,tolerance='1e-8;1e-10;1e-12',precision='float64',warmup=2,repetitions=31,order='pseudorandom'))
    rows.append(dict(block='R-order verification',tolerance='fixed cycle count',precision='mpmath',warmup='not applicable',repetitions='not applicable',order='not timed'))
    result.append((rows,['docs/PROTOCOLS.md']))
    result.append(([dict(system=f,n=n,mu0=CostModel.for_field(f).mu0,mu1_kappa1=CostModel.for_field(f).mu1,
                         mu1_expression=['2+1/n','2+30/n','2+(1+kappa)/n','2+(15+kappa)/n','2+(17+kappa)/n'][i]) for i,(f,n) in enumerate(DIMENSIONS.items())],['models.field_costs']))
    result.append(([dict(operation=a,equivalent_products=b) for a,b in [('multiplication',1),('division','kappa'),('square root',3),('exp',16),('log',14),('sin',14),('cos',14),('arctan',14)]],['Manuscript Table 6; declared weights, not new timings']))
    result.append((read_csv('structural_summary.csv'),['structural_summary.csv']))
    with (ROOT/'data/manuscript/initialization_summary.csv').open(newline='') as f: rows=list(csv.DictReader(f))
    result.append((rows,['data/manuscript/initialization_summary.csv; transcribed from Table 8, not raw samples']))
    h=next(r for r in read_csv('hammerstein_stationary_model.csv') if float(r['kappa'])==1)
    result.append(([dict(pair=f'P{N}/M{M}',cost_A=h[f'cost_P{N}'],cost_B=h[f'cost_M{M}'],eta_A=h[f'eta_P{N}'],eta_B=h[f'eta_M{M}'],advantage_A_pct=h[f'adv_eta_P{N}_pct']) for N,M in [(5,6),(7,8)]],['hammerstein_stationary_model.csv']))
    result.append((read_csv('hammerstein_pairs.csv'),['hammerstein_pairs.csv']))
    rows=[]
    for field in ('H1','H5'):
        for k in (1.0,4.45444):
            m=CostModel.for_field(field,kappa=k);N,_=m.optimal('P');s,_=m.optimal('S')
            rows.append(dict(system=field,kappa=k,N_star=N,m_star=s,eta_P=m.eta_p(N),eta_S=m.eta_s(s),
                             advantage_P_pct=100*(m.eta_p(N)/m.eta_s(s)-1),P_band=str(m.near_optimal('P')),S_band=str(m.near_optimal('S'))))
    result.append((rows,['kappa_pn_shamanskii.csv','models.CostModel.optimal']))
    result.append(([r for r in read_csv('continuous_curvature.csv') if float(r['kappa'])==1],['continuous_curvature.csv']))
    sweep=read_csv('pn_shamanskii_results.csv')
    rows=[]
    for field in ('H1','H5'):
        m=CostModel.for_field(field)
        for delta in (.1,.3,.6):
            for fam in ('P','S'):
                index=m.optimal(fam)[0]
                r=next(r for r in sweep if r['system']==field and math.isclose(float(r['delta']),delta) and r['family']==fam and int(r['index'])==index)
                rows.append({k:r[k] for k in ['system','n','delta','family','index','cycles','time_median_s','time_iqr_s']})
    result.append((rows,['pn_shamanskii_results.csv']))
    rows=[]
    for field in ('H1','H5'):
        for delta in (.1,.3,.6):
            for family in ('P','S'):
                values=[r for r in sweep if r['system']==field and math.isclose(float(r['delta']),delta) and r['family']==family]
                best=min(values,key=lambda r:float(r['time_median_s']))
                bestw=min(values,key=lambda r:float(r['work_total_k1p0']))
                rows.append(dict(system=field,delta=delta,family=family,best_by_time=best['index'],time_s=best['time_median_s'],iqr_s=best['time_iqr_s'],
                                 best_by_work=bestw['index'],work=bestw['work_total_k1p0'],first_single_cycle=min(int(r['index']) for r in values if int(r['cycles'])==1)))
    result.append((rows,['pn_shamanskii_results.csv']))
    rows=[]
    for N,M,name in [(5,6,'p5_m6_indices.csv'),(7,8,'p7_m8_indices.csv')]:
        for r in read_csv(name):
            if float(r['kappa'])==1:
                rows.append(dict(pair=f'P{N}/M{M}',system=r['system'],n=r['n'],eta_A=r[f'eta_P{N}'],eta_B=r[f'eta_M{M}'],advantage_A_pct=100*float(r[f'relative_eta_P{N}_vs_M{M}'])))
    result.append((rows,['p5_m6_indices.csv','p7_m8_indices.csv']))
    rows=[]
    for N,M,name in [(5,6,'p5_m6_pairs.csv'),(7,8,'p7_m8_pairs.csv')]:
        for field in DIMENSIONS:
            rr=[r for r in read_csv(name) if r['system']==field]
            times=[100*float(r[f'CPU_relative_P{N}_vs_M{M}']) for r in rr]
            work=[100*float(r[f'work_relative_P{N}_vs_M{M}']) for r in rr]
            rows.append(dict(pair=f'P{N}/M{M}',system=field,cases=len(rr),cpu_median_pct=statistics.median(times),cpu_min_pct=min(times),cpu_max_pct=max(times),
                             work_min_pct=min(work),work_max_pct=max(work),resolved=sum(int(r['CPU_difference_clear_2IQR']) for r in rr)))
    result.append((rows,['p5_m6_pairs.csv','p7_m8_pairs.csv']))
    result.append((read_csv('fusion_resources.csv'),['fusion_resources.csv']))
    result.append((read_csv('fusion_stationary_pairs.csv'),['fusion_stationary_pairs.csv; see KNOWN_ISSUES.md for G5 work coefficients']))
    rows=[]
    for family in ('G1','G5'):
        for n in (40,55,100,200):
            rr=[r for r in read_csv('fusion_end_to_end_pairs.csv') if r['family']==family and int(r['n'])==n]
            times=[100*float(r['CPU_relative_hat_vs_comp']) for r in rr]
            works=[100*float(r['work_relative_hat_vs_comp']) for r in rr]
            rows.append(dict(family=family,n=n,cpu_median_pct=statistics.median(times),cpu_min_pct=min(times),cpu_max_pct=max(times),
                             work_relative_pct=statistics.median(works),resolved=sum(int(r['CPU_clear_2IQR']) for r in rr)))
    result.append((rows,['fusion_end_to_end_pairs.csv; see KNOWN_ISSUES.md for G5 work coefficients']))
    result.append((read_csv('phat2_p10_field_summary.csv'),['phat2_p10_field_summary.csv']))
    result.append((read_csv('high_precision_orders.csv'),['high_precision_orders.csv']))
    result.append((read_csv('q_order_verification.csv'),['q_order_verification.csv']))
    assert len(result)==22
    return result


def escape(text):
    return str(text).replace('\\',r'\textbackslash{}').replace('_',r'\_').replace('%',r'\%').replace('&',r'\&').replace('#',r'\#')


def export_tables(output: Path) -> list[dict]:
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    manifest=[]
    for number,((label,title),(rows,sources)) in enumerate(zip(TABLES,table_rows()),1):
        path=output/f'table_{number:02d}.csv'
        write_csv(path,rows)
        # Full precision machine-readable exports, not editorially rounded tables.
        manifest.append(dict(table=number,label=label,title=title,file=path.name,sources=sources,rows=len(rows)))
    write_json(output/'table_manifest.json',manifest)
    return manifest
