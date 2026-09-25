"""Command-line entry points for verification, reporting, and fresh experiments."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import subprocess
import sys
from .io import ROOT


def main():
    parser=argparse.ArgumentParser(description='Reproducibility, validation, and benchmark tools for predictor-memory frozen-Jacobian methods.')
    sub=parser.add_subparsers(dest='command',required=True)
    verify=sub.add_parser('verify',help='Check archived data, algebraic formulas, trajectories and resource counts')
    verify.add_argument('--outdir',type=Path,default=Path('build/verification'))
    verify.add_argument('--models-only',action='store_true')
    reproduce=sub.add_parser('reproduce',help='Rebuild all English table exports, figures and the computational supplement from archived data')
    reproduce.add_argument('--outdir',type=Path,default=Path('build/reproduced'))
    reproduce.add_argument('--compile-pdf',action='store_true',help='Compile the supplement with pdflatex')
    benchmark=sub.add_parser('benchmark',help='Run a fresh timing experiment; never replace archived measurements')
    benchmark.add_argument('--suite',choices=['all','initialization','sweep','p5-m6','p7-m8','hammerstein','fusion','stationary','phat2-p10'],default='all')
    benchmark.add_argument('--outdir',type=Path,default=Path('build/current'))
    benchmark.add_argument('--repeats',type=int,default=None)
    benchmark.add_argument('--warmup',type=int,default=None)
    benchmark.add_argument('--threads',type=int,default=None,help='Explicit BLAS thread limit; default preserves the environment')
    benchmark.add_argument('--quick',action='store_true',help='Smoke test only: first group, up to 3 repetitions; not a publication benchmark')
    models=sub.add_parser('models',help='Regenerate parameter selection, curvature and interval checks')
    models.add_argument('--outdir',type=Path,default=Path('build/models'))
    structure=sub.add_parser('structure',help='Recompute ranks, condition numbers and complex-step Jacobian checks')
    structure.add_argument('--outdir',type=Path,default=Path('build/structure'))
    orders=sub.add_parser('orders',help='Run the original high-precision COC protocols in English (not timed)')
    orders.add_argument('--outdir',type=Path,default=Path('build/high_precision'))
    sub.add_parser('release-check',help='Verify every packaged file against SHA256SUMS')
    sub.add_parser('write-manifest',help='Maintainer command: record a new manifest only after reviewing intentional changes')
    args=parser.parse_args()
    if args.command=='verify':
        from .validation import run_validation,audit_archived_work
        result=run_validation(args.outdir/'validation.json',not args.models_only)
        result['archived_work_audit']=audit_archived_work(args.outdir)
    elif args.command=='reproduce':
        from .reporting import export_tables
        from .figures import create_figures
        from .supplement import generate_supplement
        tables=export_tables(args.outdir/'tables');figures=create_figures(args.outdir/'figures')
        supplement=generate_supplement(args.outdir/'supplement',args.outdir/'figures',args.compile_pdf)
        result=dict(tables=len(tables),figure_panels=len(figures),supplement=str(supplement))
    elif args.command=='benchmark':
        from .benchmarks import run_benchmark
        result=run_benchmark(args.suite,args.outdir,repeats=args.repeats,warmup=args.warmup,threads=args.threads,quick=args.quick)
    elif args.command=='models':
        from .analysis import parameter_studies
        result=parameter_studies(args.outdir)
    elif args.command=='structure':
        from .analysis import structural_audit
        result=structural_audit(args.outdir)
    elif args.command=='orders':
        for script in ('verify_orders.py','verify_q_orders.py'):
            subprocess.run([sys.executable,str(ROOT/'experiments'/script),'--outdir',str(args.outdir.resolve())],check=True)
        from .validation import check_high_precision_output
        result=check_high_precision_output(args.outdir)
    else:
        from .maintenance import verify_manifest,write_manifest
        result=write_manifest() if args.command=='write-manifest' else verify_manifest()
    print(json.dumps(result,indent=2,ensure_ascii=False))

if __name__=='__main__': main()
