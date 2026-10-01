"""Paired and resumable multiprecision experiments."""
from __future__ import annotations
import csv, gzip, hashlib, json, os, platform, random, sys, time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import mpmath
from mpmath import mp
from mpmath.libmp import BACKEND
from .core import DenseProblem, StressProblem, FIELDS, Operations, coefficients, factor, solve, norm, text, pack_vec, unpack_vec, method, iterate, traces, step, expected
from .costs import optimum, order, eta, work
ROOT=Path(__file__).resolve().parents[1]


def stamp(): return datetime.now(timezone.utc).isoformat()

def digest(obj): return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def fhash(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def code_hash():
    paths=[ROOT/'run_family_multiprecision.py']+sorted((ROOT/'pn_family_mp').glob('*.py'))
    return digest({str(p.relative_to(ROOT)):fhash(p) for p in paths})


def js(path,obj):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n',encoding='utf-8'); tmp.replace(path)


def csvout(path,rows):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8',newline='') as f:
        if rows:
            w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def gzout(path,obj):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    with gzip.open(path,'wt',encoding='utf-8') as f: json.dump(obj,f,ensure_ascii=False)


def env():
    try:
        import gmpy2
        g={'version':gmpy2.version(),'gmp':gmpy2.mp_version()}
    except ImportError: g=None
    return {'utc':stamp(),'platform':platform.platform(),'processor':platform.processor(),
            'machine':platform.machine(),'python':sys.version,'executable':sys.executable,
            'numpy':np.__version__,'mpmath':mpmath.__version__,'mpmath_backend':BACKEND,
            'gmpy2':g,'reported_logical_cpus':os.cpu_count(),'command':sys.argv,
            'timer':vars(time.get_clock_info('perf_counter')),'code_sha256':code_hash(),
            'thread_environment':{k:os.environ.get(k) for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')},
            'arithmetic':'Serial mpmath mpf / NumPy object arrays; no BLAS in the MP solver',
            'timer_scope':'Full solve including copies, counters and terminal residuals; excluding input setup, error norms, trace serialization and output I/O'}


def gid(g): return f"{g['suite']}_{g['field']}_n{g['n']}_d{g['delta'].replace('.','p')}_e{g['digits']}"


def proposed_plan():
    groups=[]
    for name in ('H1','H5'):
        n=FIELDS[name]
        methods=[f"P{optimum('P',name,n)}",f"S{optimum('S',name,n)}"]
        for t in [50,200,800,1024]:
            groups.append({'suite':'families','field':name,'n':n,'delta':'0.30','digits':t,'methods':methods})
    return {'schema':1,'status':'PROSPECTIVE',
            'working_dps':1300,'verification_dps':2600,'max_cycles':20,'warmups':2,
            'repetitions':31,'bootstrap_resamples':10000,'seed':20260925,
            'initialization':'D','reference_kappa':1,
            'cost_scope':'Canonical reference indices, not calibrated multiprecision time optima',
            'input_rule':'Frozen float64 A/B/C lifted as exact dyadics; analytic constants/alpha in MP; one working-precision x0 exactly lifted for validation',
            'gate_relative_limit':'1e-10','gate_guard_digits':40,
            'groups':groups}


def validate_plan(plan):
    groups=plan['groups']
    if not groups: raise ValueError('Empty experiment plan')
    if plan['working_dps']<max(g['digits'] for g in groups)+150: raise ValueError('At least 150 guard digits required')
    if plan['verification_dps']<2*plan['working_dps']: raise ValueError('Verification dps must be at least twice working dps')
    if plan['warmups']<1 or plan['repetitions']<7: raise ValueError('At least 1 warmup and 7 repetitions; defaults 2/31')
    if len({gid(g) for g in groups})!=len(groups): raise ValueError('Duplicate group identifiers')
    for g in groups:
        if g['field'] not in FIELDS or g['n']<2 or g['digits']<1: raise ValueError('Invalid problem or tolerance')
        if len(g['methods'])!=len(set(g['methods'])): raise ValueError('Duplicate methods')
        for key in g['methods']: method(key)
        if g['suite']=='families':
            correct=[f"P{optimum('P',g['field'],g['n'])}",f"S{optimum('S',g['field'],g['n'])}"]
            if g['methods']!=correct: raise ValueError('Family comparison must use independent reference-index optima')


def input_name(n,d): return f"x0_n{n}_d{d.replace('.','p')}.json"


def prepare(out,plan):
    directory=Path(out)/'inputs'; directory.mkdir(parents=True,exist_ok=True)
    for n in sorted({g['n'] for g in plan['groups']}):
        path=directory/f'coefficients_n{n}.npz'
        if not path.exists():
            A,B,C=coefficients(n); np.savez_compressed(path,A=A,B=B,C=C)
    identities={}
    with mp.workdps(plan['working_dps']):
        for n,d in sorted({(g['n'],g['delta']) for g in plan['groups']}):
            p=DenseProblem('H1',n,directory); identities[str(n)]=p.input_sha256
            path=directory/input_name(n,d)
            value={'dps':mp.dps,'n':n,'delta':d,'x':pack_vec(p.initial(d))}
            if path.exists() and json.loads(path.read_text())!=value: raise ValueError('An incompatible x0 already exists; use a new output directory')
            js(path,value)
    js(directory/'identity.json',identities)
    return input_hash(out)


def input_hash(out):
    return digest({p.name:fhash(p) for p in sorted((Path(out)/'inputs').glob('*')) if p.is_file()})


def x0_from(out,g):
    return unpack_vec(json.loads((Path(out)/'inputs'/input_name(g['n'],g['delta'])).read_text())['x'])


def serial(p,r):
    return {**{k:r[k] for k in ('method','status','cycles','active_outer_cycles','active_predictors','counts')},
            'residual':text(r['residual']),'error':text(norm(r['x']-p.alpha)), 'history':traces(p,r)}


def precision_gate(out,g,key,plan):
    """Empirical numerical check, not an interval-certified error bound."""
    D,V=plan['working_dps'],plan['verification_dps']; timing={}
    with mp.workdps(D):
        p=DenseProblem(g['field'],g['n'],Path(out)/'inputs'); x0=x0_from(out,g)
        start=time.perf_counter(); low=iterate(p,x0,method(key),mp.mpf(f"1e-{g['digits']}"),plan['max_cycles'],trace=True)
        timing['low_run_s']=time.perf_counter()-start; low_doc=serial(p,low)
    with mp.workdps(V):
        p=DenseProblem(g['field'],g['n'],Path(out)/'inputs'); x0=x0_from(out,g); tol=mp.mpf(f"1e-{g['digits']}")
        start=time.perf_counter(); high=iterate(p,x0,method(key),tol,plan['max_cycles'],trace=True)
        timing['high_run_s']=time.perf_counter()-start; high_doc=serial(p,high)
        rr=norm(p.F(low['x'])); ee=norm(low['x']-p.alpha)
        floor=mp.power(10,-D+plan['gate_guard_digits'])*(1+norm(p.alpha)); diffs=[]
        for k,(xl,xh) in enumerate(zip(low['history'],high['history'])):
            e=norm(xh-p.alpha); d=norm(xl-xh)/max(e,floor)
            diffs.append({'iteration':k,'relative_difference_or_floor_scaled':text(d),
                          'floor_censored':bool(e<=floor),'pass':bool(d<mp.mpf(plan['gate_relative_limit']))})
        passed=low['status']=='converged' and high['status']=='converged' and low['cycles']==high['cycles'] and rr<=tol and all(t['pass'] for t in diffs)
        return {'group':gid(g),'method':key,'passed':bool(passed),'working_dps':D,'verification_dps':V,
                'checked_residual_of_low_iterate':text(rr),'checked_error_of_low_iterate':text(ee),
                'low_cycles':low['cycles'],'high_cycles':high['cycles'], 'active_outer_cycles':low['active_outer_cycles'],
                'iterate_comparison':diffs,'feasibility_durations_not_comparative_timings':timing,
                'low':low_doc,'high':high_doc}


def pilot(out,large=False,quick=False):
    out=Path(out)
    if (out/'protocol_frozen.json').exists(): raise ValueError('Do not place a pilot inside a frozen final run')
    plan=proposed_plan()
    if large: field,n,keys='H5',200,['P9','S10']
    else: field,n,keys='H1',20,['P3','S3','P5','M6','P7','M8','Phat2q2','Phat2q3','Phat2q4','P2x2','P10']
    g={'suite':'pilot','field':field,'n':n,'delta':'0.30','digits':1024,'methods':keys}
    if quick:
        plan['working_dps']=200; plan['verification_dps']=400; g['digits']=50
    plan['groups']=[g]; plan['status']='PILOT'; prepare(out,plan)
    js(out/'pilot_plan.json',plan); js(out/'environment.json',env())
    rows=[]
    for key in keys:
        print(f"Checking {field}, n={n}, {key}, {plan['working_dps']}/{plan['verification_dps']} digits ...",flush=True)
        start=time.perf_counter()
        try:
            result=precision_gate(out,g,key,plan); gzout(out/'trajectories'/f'{key}.json.gz',result)
            row={k:result[k] for k in ('method','passed','low_cycles','high_cycles','active_outer_cycles','checked_residual_of_low_iterate','checked_error_of_low_iterate')}
            row.update(result['feasibility_durations_not_comparative_timings'])
            row['exception']=''
        except Exception as e:
            row={'method':key,'passed':False,'low_cycles':'','high_cycles':'','active_outer_cycles':'',
                 'checked_residual_of_low_iterate':'','checked_error_of_low_iterate':'',
                 'low_run_s':'','high_run_s':'','exception':repr(e)}
        row['elapsed_s']=time.perf_counter()-start; rows.append(row)
        csvout(out/'pilot_summary.csv',rows)
        print(f"  pass={row['passed']}; cycles={row['low_cycles']}/{row['high_cycles']}; elapsed={row['elapsed_s']:.2f} s",flush=True)
    js(out/'pilot_status.json',{'complete':True,'passed':all(r['passed'] for r in rows),'large':large,'quick':quick,
                              'code_sha256':code_hash(),'scope':'Numerical feasibility, not a final repeated timing comparison'})
    return rows


def freeze(out,plan_path=None,neutral_note=''):
    out=Path(out); path=out/'protocol_frozen.json'
    if path.exists(): raise ValueError('Protocol already frozen; do not overwrite')
    if (out/'sessions').exists(): raise ValueError('Cannot freeze over prior observations')
    plan=proposed_plan() if plan_path is None else json.loads(Path(plan_path).read_text())
    validate_plan(plan)
    if plan_path is not None and digest(plan['groups'])!=digest(proposed_plan()['groups']) and not neutral_note.strip():
        raise ValueError('A modified grid requires a prospectively stated neutral resource/design reason')
    plan['status']='FROZEN_BEFORE_TIMINGS'; plan['freeze_utc']=stamp()
    plan['neutral_design_note']=neutral_note or 'Frozen experiment design'
    plan['code_sha256']=code_hash(); plan['input_sha256']=prepare(out,plan)
    plan['protocol_sha256']=digest(plan)
    js(path,plan); js(out/'freeze_environment.json',env())
    return plan


def frozen(out):
    p=json.loads((Path(out)/'protocol_frozen.json').read_text())
    wanted=p.pop('protocol_sha256')
    if digest(p)!=wanted: raise ValueError('Frozen protocol was modified')
    p['protocol_sha256']=wanted
    if p['code_sha256']!=code_hash(): raise ValueError('Code changed since freeze; use a new protocol/output')
    if p['input_sha256']!=input_hash(out): raise ValueError('Frozen inputs changed')
    return p


def describe(plan):
    return {'groups':len(plan['groups']),'method_configurations':sum(len(g['methods']) for g in plan['groups']),
            'timed_full_solves':plan['repetitions']*sum(len(g['methods']) for g in plan['groups']),
            'working_dps':plan['working_dps'],'verification_dps':plan['verification_dps']}


def summary_rows(raw,plan,g):
    data={key:[] for key in g['methods']}
    for row in raw: data[row['method']].append(row['elapsed_ns']*1e-9)
    count=len(next(iter(data.values())))
    if any(len(t)!=count for t in data.values()) or count!=plan['repetitions']: raise ValueError('Incomplete paired group')
    ref=g['methods'][-1]; b=np.asarray(data[ref]); rows=[]
    rng=np.random.default_rng(int(digest([plan['seed'],gid(g),'bootstrap'])[:16],16))
    inds=rng.integers(0,count,(plan['bootstrap_resamples'],count))
    for key in g['methods']:
        a=np.asarray(data[key]); ma=float(np.median(a)); mb=float(np.median(b))
        ia=float(np.percentile(a,75)-np.percentile(a,25)); ib=float(np.percentile(b,75)-np.percentile(b,25))
        boot=np.median(a[inds],axis=1)/np.median(b[inds],axis=1)
        low,high=np.percentile(boot,[2.5,97.5]); example=next(r for r in raw if r['method']==key)
        c=example['counts']; rows.append({'group':gid(g),'suite':g['suite'],'field':g['field'],'n':g['n'],
            'delta':g['delta'],'tolerance_digits':g['digits'],'method':key,'baseline':ref,'repetitions':count,
            'median_s':ma,'iqr_s':ia,'median_ratio':ma/mb,'paired_bootstrap_ratio_low':float(low),
            'paired_bootstrap_ratio_high':float(high),'two_iqr_filter':bool(abs(ma-mb)>2*(ia+ib)),
            'cycles':example['cycles'],'active_outer_cycles':example['active_outer_cycles'],
            'reference_eta':eta(key,g['field'],g['n']),'reference_work':work(c,g['field'],g['n']),**c})
    return rows


def run(out,suite=None,max_groups=None):
    out=Path(out); plan=frozen(out)
    session=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    sd=out/'sessions'/session; js(sd/'environment.json',env()); done=out/'completed'; done.mkdir(exist_ok=True)
    groups=[g for g in plan['groups'] if suite is None or g['suite']==suite]
    completed=0
    for g in groups:
        id=gid(g)
        if (done/f'{id}.json').exists(): continue
        gd=sd/id; gd.mkdir(parents=True)
        print(f'Precision gate: {id}',flush=True)
        gates={}; failed=False
        for key in g['methods']:
            try:
                r=precision_gate(out,g,key,plan); gates[key]=r; gzout(gd/f'gate_{key}.json.gz',r)
                if not r['passed']: failed=True
            except Exception as e:
                failed=True; js(gd/f'gate_{key}_error.json',{'error':repr(e)})
        if failed:
            js(gd/'status.json',{'status':'GATE_FAILED_NO_FINAL_TIMINGS','group':g,'session':session})
            # A failed group remains in the record; no automatic precision adaptation or silent deletion.
            print('Gate failed; examine recorded data before any prospective protocol amendment.',flush=True)
            continue
        with mp.workdps(plan['working_dps']):
            p=DenseProblem(g['field'],g['n'],out/'inputs'); x0=x0_from(out,g); tol=mp.mpf(f"1e-{g['digits']}")
            for _ in range(plan['warmups']):
                for key in g['methods']: iterate(p,x0,method(key),tol,plan['max_cycles'])
            raw=[]; rng=random.Random(int(digest([plan['seed'],id])[:16],16)); logfile=gd/'observations.jsonl'
            for rep in range(plan['repetitions']):
                keys=list(g['methods']); rng.shuffle(keys)
                for position,key in enumerate(keys):
                    start=time.perf_counter_ns(); r=iterate(p,x0,method(key),tol,plan['max_cycles']); elapsed=time.perf_counter_ns()-start
                    valid=r['status']=='converged' and r['cycles']==gates[key]['low_cycles'] and r['counts']==gates[key]['low']['counts']
                    row={'session':session,'group':id,'repeat':rep,'position':position,'method':key,
                         'elapsed_ns':elapsed,'valid':valid,'status':r['status'],'cycles':r['cycles'],
                         'active_outer_cycles':r['active_outer_cycles'],'active_predictors':r['active_predictors'],
                         'counts':r['counts'],'residual':text(r['residual']),'error':text(norm(r['x']-p.alpha))}
                    with logfile.open('a',encoding='utf-8') as f: f.write(json.dumps(row)+'\n'); f.flush()
                    raw.append(row)
                    if not valid: raise ArithmeticError('Timing result diverged from precision-gated trajectory; group is incomplete')
                print(f'  {id}: repetition {rep+1}/{plan["repetitions"]}',flush=True)
            rows=summary_rows(raw,plan,g); csvout(gd/'summary.csv',rows)
            js(done/f'{id}.json',{'status':'COMPLETE','session':session,'directory':str(gd.relative_to(out)),
                                  'protocol_sha256':plan['protocol_sha256'],'group':g})
        completed+=1
        if max_groups is not None and completed>=max_groups: break
    summarize(out)


def summarize(out):
    out=Path(out); plan=frozen(out); rows=[]
    for p in sorted((out/'completed').glob('*.json')):
        d=json.loads(p.read_text())
        if d['protocol_sha256']!=plan['protocol_sha256']: raise ValueError('Mixed protocol hashes')
        with (out/d['directory']/'summary.csv').open(newline='',encoding='utf-8') as f:
            for row in csv.DictReader(f): row['session']=d['session']; rows.append(row)
    csvout(out/'combined_summary.csv',rows)
    completed=len(list((out/'completed').glob('*.json')))
    js(out/'coverage.json',{'planned_groups':len(plan['groups']),'complete_groups':completed,
                          'all_groups_complete':completed==len(plan['groups']),
                          'incomplete_policy':'Only complete within-session paired groups are aggregated; partial/failed groups remain under sessions'})
    return rows
