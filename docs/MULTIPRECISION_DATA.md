# Multiprecision datasets

This repository contains two high-precision timing datasets for the frozen-Jacobian study.

## Family comparison

The family dataset compares the independently selected reference-index optima of the predictor-memory family and Shamanskii's family on two dense test fields:

- `H1`, dimension 20: `P3` versus `S3`;
- `H5`, dimension 200: `P9` versus `S10`.

The starting distance is 0.30, the residual tolerances are `1e-50`, `1e-200`, `1e-800`, and `1e-1024`, the working precision is 1300 decimal digits, and the verification precision is 2600 digits. There are 8 groups, 496 timed observations, and 16 higher-precision controls.

## External comparators

The external dataset compares `P5/M6` and `P7/M8` on `H1` and `H5` using the same starting distance, tolerances, and precision levels. Each group begins with 15 paired repetitions. A predeclared mechanical rule extends an ambiguous group by exactly 16 additional pairs, giving 31. There are 16 final groups, 736 timed observations, and 48 higher-precision controls.

Platform records are preserved with the data. AMD and Intel observations are not pooled within a final comparison group.

## Interpretation

The repository keeps separate:

1. local/asymptotic convergence order;
2. the algebraic reference-efficiency index;
3. work required to reach a prescribed tolerance; and
4. measured elapsed wall time for a specified implementation and platform.

The timing data do not establish a machine-independent performance ordering. The datasets compare complete algorithms and do not isolate a causal effect attributable only to memory.

The multiprecision timing datasets in this repository do not include a timing-sensitivity study of the fused family for `q=2,3,4`.
