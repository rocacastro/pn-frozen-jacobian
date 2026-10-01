#!/usr/bin/env python3
"""English command-line interface. Run --help before a final timing campaign."""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
# mpmath may internally convert large mantissas when printing tiny diagnostics.
# The source and inputs here are local; use a finite, explicit bound.
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(max(sys.get_int_max_str_digits(), 1000000))
from pn_family_mp import experiments as ex


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    sub=ap.add_subparsers(dest='command',required=True)
    p=sub.add_parser('plan'); p.add_argument('--save',default='configs/protocol.json')
    p=sub.add_parser('pilot'); p.add_argument('--out',default='results/pilot_small'); p.add_argument('--large',action='store_true'); p.add_argument('--quick',action='store_true')
    p=sub.add_parser('freeze'); p.add_argument('--out',default='results/final'); p.add_argument('--plan'); p.add_argument('--reason',default='')
    p=sub.add_parser('run'); p.add_argument('--out',default='results/final'); p.add_argument('--suite',choices=['families','external6','external8','fusion','secondary']); p.add_argument('--max-groups',type=int)
    p=sub.add_parser('summarize'); p.add_argument('--out',default='results/final')
    p=sub.add_parser('test'); p.add_argument('--out',default='validation/self_test.json')
    p=sub.add_parser('order'); p.add_argument('--out',default='results/order'); p.add_argument('--quick',action='store_true')
    p=sub.add_parser('kernels'); p.add_argument('--out',default='results/kernels'); p.add_argument('--field',default='H1',choices=list(ex.FIELDS)); p.add_argument('--n',type=int); p.add_argument('--dps',type=int,default=1300); p.add_argument('--repeats',type=int,default=31)
    p=sub.add_parser('stationary'); p.add_argument('--out',default='results/stationary'); p.add_argument('--field',default='H1',choices=['H1','H5']); p.add_argument('--dps',type=int,default=1300); p.add_argument('--repeats',type=int,default=31)
    args=ap.parse_args()
    if args.command=='plan':
        p=ex.proposed_plan(); ex.js(args.save,p); print(json.dumps(ex.describe(p),indent=2)); print(f'Saved experiment plan: {args.save}')
    elif args.command=='pilot': ex.pilot(args.out,args.large,args.quick)
    elif args.command=='freeze': print(json.dumps(ex.describe(ex.freeze(args.out,args.plan,args.reason)),indent=2))
    elif args.command=='run':
        if args.max_groups is not None and args.max_groups<1: raise ValueError('--max-groups must be positive')
        ex.run(args.out,args.suite,args.max_groups)
    elif args.command=='summarize': ex.summarize(args.out)
    elif args.command=='test':
        from pn_family_mp.validation import self_test
        self_test(args.out)
    else:
        from pn_family_mp import diagnostics
        if args.command=='order': diagnostics.order_checks(args.out,args.quick)
        elif args.command=='kernels': diagnostics.kernels(args.out,args.field,args.n,args.dps,args.repeats)
        elif args.command=='stationary': diagnostics.stationary(args.out,args.field,args.dps,args.repeats)

if __name__=='__main__':
    try: main()
    except KeyboardInterrupt:
        print('\nInterrupted. Partial observations are preserved. Resume with the same command; only complete groups are skipped.',file=sys.stderr)
        raise SystemExit(130)
    except Exception as e:
        print(f'ERROR: {type(e).__name__}: {e}',file=sys.stderr)
        raise
