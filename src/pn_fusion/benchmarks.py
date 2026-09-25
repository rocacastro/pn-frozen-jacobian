"""Fresh, interleaved timing experiments; archived publication timings are immutable."""
from __future__ import annotations
from dataclasses import dataclass,asdict
from datetime import datetime, timezone
from pathlib import Path
import contextlib
import hashlib
import math
import os
import platform
import statistics
import time
import numpy as np
import scipy
from threadpoolctl import threadpool_info,threadpool_limits
from .io import write_csv,write_json
from .models import DIMENSIONS,CostModel,rho,fused_order
from .problems import dense_problem,hammerstein_problem
from .methods import solve_problem,prepare_stationary,stationary_macrocycle

DELTAS=(0.1,0.3,0.6)
TOLERANCES=(1e-8,1e-10,1e-12)

@dataclass(frozen=True)
class Configuration:
    name: str
    method: str
    N: int=2
    q: int=4
    initialization: str='D'


def capture_environment() -> dict:
    libraries=[]
    for info in threadpool_info():
        libraries.append({key:value for key,value in info.items() if key!='filepath'})
    return dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),platform=platform.platform(),
                processor=platform.processor(),python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                libraries=libraries,clock='time.perf_counter_ns; elapsed wall-clock time',
                timing_scope='Fresh run of the English reimplementation, not a replacement for archived timings',
                initialization='D except for explicit initialization comparison',maxit=20)


def groups(suite: str):
    """Yield predeclared groups; no parameter depends on observed timings."""
    if suite=='sweep':
        for field in ('H1','H5'):
            for delta in DELTAS:
                configs=[Configuration(f'P{N}','P',N) for N in range(2,21)]
                configs += [Configuration(f'S{m}','S',m) for m in range(1,21)]
                yield field,DIMENSIONS[field],delta,1e-12,configs,False
    elif suite in ('p5-m6','p7-m8'):
        N,external=(5,'M6') if suite=='p5-m6' else (7,'M8')
        for field,n in DIMENSIONS.items():
            for delta in DELTAS:
                for tol in TOLERANCES:
                    yield field,n,delta,tol,[Configuration(f'P{N}','P',N),Configuration(external,external)],False
    elif suite=='initialization':
        cases=[('H1',3),('H5',9)]+[(field,N) for N in (5,7) for field in DIMENSIONS]
        for field,N in cases:
            yield field,DIMENSIONS[field],0.1,1e-12,[Configuration('I','P',N,initialization='I'),Configuration('D','P',N)],False
    elif suite=='hammerstein':
        for tol in TOLERANCES:
            yield 'Hammerstein',8,0.0,tol,[Configuration('P5','P',5),Configuration('M6','M6'),
                                         Configuration('P7','P',7),Configuration('M8','M8')],False
    elif suite=='phat2-p10':
        for field in ('H3','H5'):
            for delta in DELTAS:
                for tol in TOLERANCES:
                    yield field,DIMENSIONS[field],delta,tol,[Configuration('Phat2','Phat'),Configuration('P10','P',10)],False
    elif suite in ('fusion','stationary'):
        for field in ('G1','G5'):
            for n in (40,55,100,200):
                for delta in ((0.1,) if suite=='stationary' else DELTAS):
                    for tol in ((1e-12,) if suite=='stationary' else TOLERANCES):
                        yield field,n,delta,tol,[Configuration('Phat2','Phat'),Configuration('P2composeP2','Pcompose')],suite=='stationary'
    else:
        raise ValueError(f'Unknown benchmark suite: {suite}')


def run_benchmark(suite: str, output: Path, *, repeats: int|None=None,
                  warmup: int|None=None, threads: int|None=None, quick: bool=False) -> dict:
    suites=['initialization','sweep','p5-m6','p7-m8','hammerstein','fusion','stationary','phat2-p10'] if suite=='all' else [suite]
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    if threads is not None and threads<1: raise ValueError('threads must be positive')
    all_reports=[]
    with (threadpool_limits(limits=threads,user_api='blas') if threads is not None else contextlib.nullcontext()):
        environment=capture_environment()
        environment.update(requested_threads=threads,quick_mode=quick)
        write_json(output/'environment.json',environment)
        for current in suites:
            reps=repeats if repeats is not None else (51 if current=='initialization' else 31)
            warm=warmup if warmup is not None else (3 if current=='initialization' else 2)
            if reps<1 or warm<0: raise ValueError('repeats must be >= 1 and warmup >= 0')
            if quick: reps,warm=min(reps,3),min(warm,1)
            raw_rows=[]; result_rows=[]; pair_rows=[]
            all_groups=list(groups(current))
            if quick: all_groups=all_groups[:1]
            for group_index,(field,n,delta,tol,configs,stationary) in enumerate(all_groups):
                problem=hammerstein_problem() if field=='Hammerstein' else dense_problem(field,delta,n)
                model=CostModel.for_field(field,n)
                states={c.name:prepare_stationary(problem,c.method,c.N,c.q) for c in configs} if stationary else {}
                def execute(c):
                    if stationary: return stationary_macrocycle(problem,c.method,states[c.name],c.N,c.q)
                    return solve_problem(problem,c.method,N=c.N,q=c.q,tol=tol,initialization=c.initialization)
                for _ in range(warm):
                    for c in configs: execute(c)
                seed=20260923+n+int(round(delta*1000))+int(round(-math.log10(tol)))
                rng=np.random.Generator(np.random.PCG64(seed))
                times={c.name:[] for c in configs}; last={}; signatures={}
                for repetition in range(reps):
                    order=list(configs)
                    if current=='initialization':
                        if repetition%2: order.reverse()
                    else: rng.shuffle(order)
                    for position,c in enumerate(order):
                        start=time.perf_counter_ns();answer=execute(c);elapsed=(time.perf_counter_ns()-start)*1e-9
                        signature=(answer.cycles,answer.converged,tuple(asdict(answer.resources).values()))
                        if c.name in signatures and signature!=signatures[c.name]:
                            raise RuntimeError('Nondeterministic trajectory/resource counts within one timing group')
                        signatures[c.name]=signature;times[c.name].append(elapsed);last[c.name]=answer
                        raw_rows.append(dict(suite=current,system=field,n=n,delta=delta,tol=tol,method=c.name,
                                             repetition=repetition,position=position,elapsed_s=elapsed,shuffle_seed=seed))
                summaries={}
                for c in configs:
                    ans=last[c.name];v=times[c.name]
                    median=statistics.median(v);iqr=float(np.percentile(v,75)-np.percentile(v,25))
                    row=dict(suite=current,system=field,n=n,delta=delta,tol=tol,method=c.name,N=c.N,q=c.q,
                             initialization=c.initialization,stationary_only=int(stationary),cycles=ans.cycles,
                             converged=int(ans.converged),residual=ans.residual,error=ans.error,
                             F_total_including_terminal=ans.resources.F,**{k:value for k,value in asdict(ans.resources).items() if k!='F'},
                             work_kappa1=ans.resources.work(model),time_median_s=median,time_iqr_s=iqr,
                             time_min_s=min(v),time_max_s=max(v),repeats=reps,warmup=warm)
                    result_rows.append(row);summaries[c.name]=row
                pairs=[(configs[0].name,configs[1].name)] if len(configs)==2 else []
                if current=='hammerstein': pairs=[('P5','M6'),('P7','M8')]
                for A,B in pairs:
                    a,b=summaries[A],summaries[B]
                    pair_rows.append(dict(suite=current,system=field,n=n,delta=delta,tol=tol,pair=A+'/'+B,
                        cycles_A=a['cycles'],cycles_B=b['cycles'],both_converged=int(a['converged'] and b['converged']),cpu_relative_A_pct=100*(a['time_median_s']/b['time_median_s']-1),
                        work_relative_A_pct=100*(a['work_kappa1']/b['work_kappa1']-1),
                        clear_2iqr=int(a['converged'] and b['converged'] and abs(a['time_median_s']-b['time_median_s'])>2*(a['time_iqr_s']+b['time_iqr_s']))))
                # Save after every group so a long run can be inspected safely.
                write_csv(output/f'{current}_samples.csv',raw_rows)
                write_csv(output/f'{current}_results.csv',result_rows)
                if pair_rows: write_csv(output/f'{current}_pairs.csv',pair_rows)
                print(f'{current}: group {group_index+1}/{len(all_groups)} completed',flush=True)
            all_reports.append(dict(suite=current,configurations=len(result_rows),raw_samples=len(raw_rows),warmup=warm,repeats=reps))
        write_json(output/'run_summary.json',all_reports)
    return {'runs':all_reports}
