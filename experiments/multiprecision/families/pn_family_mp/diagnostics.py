"""Separate order, kernel and trajectory-prepared stationary diagnostics."""
from __future__ import annotations
import json, random, time
from pathlib import Path
import numpy as np
from mpmath import mp
from .core import *
from .costs import optimum, order
from .experiments import js, csvout, gzout, env, digest, code_hash, serial, prepare, x0_from, proposed_plan, summary_rows


def new_run(out,config):
    out=Path(out)
    if (out/'config.json').exists(): raise ValueError('Diagnostic output already exists; use a new directory')
    config={**config,'code_sha256':code_hash(),'frozen_before_timings':True}
    js(out/'config.json',config); js(out/'environment.json',env())
    return out


def order_checks(out,quick=False):
    cases=[(f'P{N}',5 if N<=6 else 4,6000) for N in range(2,9)]
    cases += [(f'Phat{N}q{N+2}',4 if N==2 else 3,30000) for N in (2,3,4)]
    cases += [('Phat2q2',4,6000),('Phat2q3',4,30000)]
    if quick: cases=[('P2',3,300),('Phat2q3',2,500)]
    out=new_run(out,{'type':'fixed_cycle_order_diagnostics','quick':quick,'cases':cases,
                      'x0':'Decimal strings constructed at working precision, then exactly lifted for validation',
                      'validation':'Twice working precision; relative trajectory discrepancy < 1e-8; zeros are not certified'})
    rows=[]
    for key,cycles,dps in cases:
        print(f'Order {key}: {cycles} cycles at {dps}/{2*dps} digits',flush=True)
        with mp.workdps(dps):
            p=StressProblem(); x0=p.initial(); low=iterate(p,x0,method(key),fixed=cycles,trace=True)
            doc=serial(p,low)
        with mp.workdps(2*dps):
            p=StressProblem(); high=iterate(p,x0,method(key),fixed=cycles,trace=True); hi=serial(p,high)
            checks=[]
            for k,(a,b) in enumerate(zip(low['history'],high['history'])):
                e=norm(b); diff=norm(a-b)
                checks.append({'iteration':k,'relative_difference':None if not e else text(diff/e),
                               'passed':bool(e and diff/e<mp.mpf('1e-8'))})
            ok=low['cycles']==cycles and high['cycles']==cycles and all(c['passed'] for c in checks)
            useful=[v['coc'] for v in doc['history'] if v['stationary_coc']]
            result={'method':key,'cycles_requested':cycles,'working_dps':dps,'verification_dps':2*dps,
                    'precision_stable':bool(ok),'comparison':checks,'low':doc,'high':hi}
            gzout(out/f'{key}_trajectory.json.gz',result)
            rows.append({'method':key,'cycles':low['cycles'],'working_dps':dps,'verification_dps':2*dps,
                         'reference_order':order(key),'last_coc':useful[-1] if useful else '',
                         'useful_stationary_coc_count':len(useful),'precision_stable':bool(ok),
                         'minus_log10_final_error':text(-mp.log10(norm(low['x']))) if norm(low['x']) else 'computed_zero'})
        csvout(out/'order_summary.csv',rows)
    return rows


def kernels(out,field='H1',n=None,dps=1300,repeats=31):
    if repeats<7 or dps<30: raise ValueError('Use at least 7 repetitions and 30 digits')
    n=FIELDS[field] if n is None else n
    out=new_run(out,{'type':'kernel_measurements','field':field,'n':n,'dps':dps,'warmups':2,'repeats':repeats,
                    'scope':'Absolute kernel wall times; empty call is a timer/call baseline, not estimated global overhead'})
    with mp.workdps(dps):
        p=DenseProblem(field,n); x=p.initial('0.30'); A=p.J(x); fac=factor(A); b=p.F(x)
        functions={'F':lambda:p.F(x),'J':lambda:p.J(x),'LU':lambda:factor(A),
                   'solve':lambda:solve(fac,b),'matvec':lambda:A@b,'empty_call':lambda:None}
        for _ in range(2):
            for fun in functions.values(): fun()
        raw=[]; rng=random.Random(9172026)
        for rep in range(repeats):
            names=list(functions); rng.shuffle(names)
            for position,key in enumerate(names):
                start=time.perf_counter_ns(); result=functions[key](); elapsed=time.perf_counter_ns()-start
                raw.append({'repeat':rep,'position':position,'kernel':key,'elapsed_ns':elapsed})
            csvout(out/'kernel_raw.csv',raw)
        rows=[]
        for key in functions:
            a=np.array([r['elapsed_ns']*1e-9 for r in raw if r['kernel']==key])
            rows.append({'field':field,'n':n,'dps':dps,'kernel':key,'repetitions':repeats,
                         'median_s':float(np.median(a)),'iqr_s':float(np.percentile(a,75)-np.percentile(a,25))})
        csvout(out/'kernel_summary.csv',rows)
    return rows


def stationary(out,field='H1',dps=1300,repeats=31):
    if repeats<7: raise ValueError('At least 7 repetitions required')
    n=FIELDS[field]; keys=[f"P{optimum('P',field,n)}",f"S{optimum('S',field,n)}"]
    out=new_run(out,{'type':'stationary_independent_optima','field':field,'n':n,'delta':'0.30',
                    'working_dps':dps,'verification_dps':2*dps,'warmups':2,'repeats':repeats,'methods':keys,
                    'preparation':'One actual delayed-start block per method; states are not forced equal',
                    'timer_scope':'One next stationary cycle, cached starting F, terminal residual test and counts object construction inside timer',
                    'interpretation':'Wall time and logarithmic error reduction, not universal temporal efficiency'})
    snapshots={}; gates=[]
    for key in keys:
        with mp.workdps(dps):
            p=DenseProblem(field,n); x0=p.initial('0.30')
            pre=iterate(p,x0,method(key),fixed=1,retain_state=True)
            low=iterate(p,x0,method(key),fixed=2,trace=True)
            snapshots[key]=(p,pre)
        with mp.workdps(2*dps):
            ph=DenseProblem(field,n); high=iterate(ph,x0,method(key),fixed=2,trace=True)
            e=norm(high['x']-ph.alpha); floor=mp.power(10,-dps+40)*(1+norm(ph.alpha))
            rel=norm(low['x']-high['x'])/max(e,floor)
            gate={'method':key,'passed':bool(low['cycles']==2 and high['cycles']==2 and rel<mp.mpf('1e-10')),
                  'relative_discrepancy_or_floor_scaled':text(rel),'low':serial(ph,low),'high':serial(ph,high)}
            gzout(out/f'gate_{key}.json.gz',gate); gates.append(gate)
    if not all(g['passed'] for g in gates): raise ArithmeticError('Stationary precision gate failed; no timings taken')
    raw=[]; reductions=[]
    with mp.workdps(dps):
        def execute(key):
            p,pre=snapshots[key]; op=Operations(p)
            x,mem=step(op,pre['x'],pre['Fx'],method(key),pre['state'])
            res=norm(op.F(x))
            return x,res,op.c.dict()
        for _ in range(2):
            for key in keys: execute(key)
        rng=random.Random(8112026)
        for rep in range(repeats):
            names=list(keys); rng.shuffle(names)
            for pos,key in enumerate(names):
                start=time.perf_counter_ns(); x,r,c=execute(key); dt=time.perf_counter_ns()-start
                raw.append({'repeat':rep,'position':pos,'method':key,'elapsed_ns':dt})
            csvout(out/'stationary_raw.csv',raw)
        for key in keys:
            p,pre=snapshots[key]; x,r,c=execute(key)
            eb,ea=norm(pre['x']-p.alpha),norm(x-p.alpha)
            times=np.array([v['elapsed_ns']*1e-9 for v in raw if v['method']==key])
            reductions.append({'method':key,'field':field,'n':n,'dps':dps,'repetitions':repeats,
                'median_s':float(np.median(times)),'iqr_s':float(np.percentile(times,75)-np.percentile(times,25)),
                'error_before':text(eb),'error_after':text(ea),'residual_before':text(pre['residual']),
                'residual_after':text(r),'decimal_error_gain':text(mp.log10(eb/ea)) if ea else 'floor_censored',**c})
        csvout(out/'stationary_summary.csv',reductions)
    return reductions
