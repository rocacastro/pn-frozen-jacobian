#!/usr/bin/env python3
"""Adaptive multiprecision timing runner for P5/M6 and P7/M8."""
from __future__ import annotations
import argparse
import json
import sys

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(max(sys.get_int_max_str_digits(), 1000000))

from pn_ext_adaptive import experiments as ex


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="command", required=True)

    p = sub.add_parser("plan")
    p.add_argument("--save", default="configs/protocol_external_adaptive.json")

    p = sub.add_parser("freeze")
    p.add_argument("--out", default="results/external_adaptive")
    p.add_argument("--plan", default="configs/protocol_external_adaptive.json")
    p.add_argument("--reason", default="")

    p = sub.add_parser("initial")
    p.add_argument("--out", default="results/external_adaptive")
    p.add_argument("--suite", choices=["external6", "external8"])
    p.add_argument("--max-groups", type=int)

    p = sub.add_parser("decide")
    p.add_argument("--out", default="results/external_adaptive")

    p = sub.add_parser("extend")
    p.add_argument("--out", default="results/external_adaptive")
    p.add_argument("--suite", choices=["external6", "external8"])
    p.add_argument("--max-groups", type=int)

    p = sub.add_parser("summarize")
    p.add_argument("--out", default="results/external_adaptive")

    p = sub.add_parser("status")
    p.add_argument("--out", default="results/external_adaptive")

    p = sub.add_parser("test")
    p.add_argument("--out", default="validation/self_test.json")

    args = ap.parse_args()
    if args.command == "plan":
        plan = ex.experiment_plan()
        ex.js(args.save, plan)
        print(json.dumps(ex.describe(plan), indent=2))
        print(f"Saved experiment plan: {args.save}")
    elif args.command == "freeze":
        print(json.dumps(ex.describe(ex.freeze(args.out, args.plan, args.reason)), indent=2))
    elif args.command == "initial":
        if args.max_groups is not None and args.max_groups < 1:
            raise ValueError("--max-groups must be positive")
        ex.run_initial(args.out, args.max_groups, args.suite)
    elif args.command == "decide":
        doc = ex.decide_extensions(args.out)
        for d in doc["decisions"]:
            target = doc['extension_target_repetitions']
            initial = doc['initial_repetitions']
            print(f"{d['group']}: {('EXTEND TO '+str(target)) if d['extend_to_31'] else ('STOP AT '+str(initial))}; {', '.join(d['decision_reasons'])}")
    elif args.command == "extend":
        if args.max_groups is not None and args.max_groups < 1:
            raise ValueError("--max-groups must be positive")
        ex.run_extensions(args.out, args.max_groups, args.suite)
    elif args.command == "summarize":
        ex.summarize(args.out)
    elif args.command == "status":
        ex.status(args.out)
    elif args.command == "test":
        ex.self_test(args.out)
        print("Self-test passed.")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrupted. Completed groups remain preserved. Resume with the same command.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as e:
        print(f"ERROR: {type(e).__name__}: {e}", file=sys.stderr)
        raise
