# Reproducibility of the multiprecision experiments

## Arithmetic and timing

The multiprecision solvers use `mpmath` numbers in NumPy object arrays, partial-pivoting LU factorization, and reused triangular solves. The solver is serial and does not use BLAS for the multiprecision linear algebra. The monotonic `perf_counter_ns` clock measures the complete solve call, including copies, counters, numerical operations, and the terminal residual test; problem construction, warmups, trace serialization, aggregation, and file output are outside the timed interval.

The family campaign uses 31 paired repetitions per group. The external campaign uses 15 paired repetitions initially and a predeclared 15-to-31 extension rule. The rule changes only the number of paired repetitions; fields, methods, starting distance, tolerances, precision, and stopping criteria are fixed independently of the observed winner.

## Precision controls

All final groups include higher-precision controls. Working and verification precision are 1300 and 2600 decimal digits, respectively. Final iterates are reevaluated at the verification precision. Precision-floor-censored values are not treated as reliable relative-error measurements for tail COC estimates.

## Data organization

- `data/multiprecision/families/` contains the family campaign.
- `data/multiprecision/external/current/` contains the final external-comparator campaign.
- `data/multiprecision/external/platform_sensitivity/` retains separate Intel follow-up data that are not pooled into the final aggregates.
- `data/multiprecision/*/audit/` contains observation ledgers, precision controls, path summaries, and descriptive diagnostics.
- `metadata/multiprecision/manifest.sha256` provides file-integrity hashes for the reader-facing update.

The Zenodo DOI `10.5281/zenodo.22950468` identifies the archived software release `v1.0.0`. The multiprecision datasets added later are not part of that archived version.
