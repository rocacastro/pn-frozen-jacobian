# Adaptive multiprecision external-comparator campaign

This package runs the prospective high-precision timing study for the external pairs
`P5/M6` and `P7/M8` on the canonical dense fields `H1` and `H5`.

The design is adaptive only in the number of paired timing repetitions. It is **not**
adaptive in fields, methods, starting distance, tolerances, precision, or stopping rule.
Those are frozen before timing.

## Frozen scientific design

- Fields: `H1` (`n=20`) and `H5` (`n=200`).
- Pairs: `P5/M6` and `P7/M8`.
- Delayed initialization `D` for `P5` and `P7`.
- Initial distance: `delta=0.30`.
- Residual tolerances: `1e-50`, `1e-200`, `1e-800`, `1e-1024`.
- Working precision: 1300 decimal digits.
- Verification precision: 2600 decimal digits.
- Two warmups per group.
- 10,000 paired percentile-bootstrap resamples.
- Every group starts with 15 paired repetitions.

## Prospectively frozen extension rule

After all 16 groups have completed 15 paired repetitions, the program mechanically
classifies each group. A group stops at 15 **only if both** conditions hold and point in
the same direction:

1. the descriptive two-IQR separation rule is met; and
2. the paired percentile-bootstrap interval for the ratio of medians excludes 1.

Otherwise that group is extended by exactly 16 additional paired repetitions, for a
total of 31. There is no manual selection after seeing the winner.

Thus the initial stage uses 480 timed full solves. If every group is ambiguous, the
maximum total is 992 timed full solves, exactly the cost of using 31 repetitions for
every group from the outset.

## Installation and validation

```console
python -m pip install -r requirements.txt
python run_external_adaptive.py test
python run_external_adaptive.py plan --save configs/protocol_external_adaptive.json
```

The `plan` command should report 16 groups, 32 method configurations, 480 initial
timed full solves, and 992 as the worst-case maximum.

## Freeze once

Use a new output directory and freeze only once:

```console
python run_external_adaptive.py freeze --plan configs/protocol_external_adaptive.json --out results/external_adaptive --reason "Prospective external-comparator tranche: H1/H5, delta 0.30, all four high-precision tolerances, 15 paired repetitions with a predeclared mechanical extension rule to 31."
```

Do not edit source files, the frozen protocol, or frozen inputs after this command.

## Initial stage: 15 paired repetitions

The program skips groups already completed, so it can be run in manageable batches.
For one group at a time:

```console
python run_external_adaptive.py initial --out results/external_adaptive --max-groups 1
```

Or, for example, complete all eight `H1` groups in one invocation by repeating the
command or using a larger `--max-groups` value.

Check progress with:

```console
python run_external_adaptive.py status --out results/external_adaptive
```

## Freeze the extension decisions

Only after all 16 initial groups are complete:

```console
python run_external_adaptive.py decide --out results/external_adaptive
```

This creates `extension_decisions.json` and `extension_decisions.csv`. The decision is
computed from the frozen rule; it is not chosen by the user.

## Extension stage

Run only the groups mechanically flagged as ambiguous. One group at a time:

```console
python run_external_adaptive.py extend --out results/external_adaptive --max-groups 1
```

Repeat until `status` reports `final_complete: true`. Groups that were clear after 15
are never rerun. Flagged groups receive 16 new paired repetitions and finish with 31.

## Final outputs

```console
python run_external_adaptive.py summarize --out results/external_adaptive
python run_external_adaptive.py status --out results/external_adaptive
```

Important files include:

- `protocol_frozen.json`
- `freeze_environment.json`
- `combined_initial_summary.csv`
- `extension_decisions.json`
- `extension_decisions.csv`
- `combined_final_summary.csv`
- `group_status.csv`
- `coverage.json`
- `initial_completed/`
- `extension_completed/`
- `sessions_initial/`
- `sessions_extension/`
- `inputs/`

Preserve the entire output directory. Pilot samples or incomplete sessions are not
silently merged into completed groups.
