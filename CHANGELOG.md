# Changelog

## [Unreleased]

### Added
- Multiprecision predictor-memory/Shamanskii family data on `H1` and `H5`: 8 groups, 496 timed observations, and 16 higher-precision controls.
- Adaptive multiprecision `P5/M6` and `P7/M8` data on `H1` and `H5`: 16 final groups, 736 timed observations, and 48 higher-precision controls.
- Frozen protocols, individual observations, phase records, platform/session metadata, precision controls, and reproducibility documentation.
- Updated computational supplement with elapsed wall-time terminology and explicit conditioning information for the dense test fields.

### Clarified
- Measurements obtained with `perf_counter_ns` are described as elapsed wall time rather than CPU time; numerical timing values are unchanged.
- The dense test Jacobians are dense and coupled but well conditioned, with observed `kappa_2(J)` approximately in `[1.404, 1.497]`.

### Archival boundary
- The tag `v1.0.0` and DOI `10.5281/zenodo.22950468` identify the initial archived software release and are not modified by this update.
- Multiprecision data added after that release are not attributed to its version-specific DOI.


## 1.0.0 — Initial public release

- Added English implementations of the iterative methods and benchmark runners used in the computational study.
- Preserved archived numerical values while standardizing dataset names and repository-facing text in English.
- Added deterministic checks for trajectories, operation counts, matrix fingerprints, and computational-cost formulas.
- Reproduced the high-precision order-verification protocols with the English implementations.
- Added reproducible exports for the article tables and the main and supplementary figures.
- Added reader-facing supplementary computational material in English.
- Added per-repetition timing capture, environment records, and repository-integrity checks for new runs.
- Documented the archived G5 evaluation-cost discrepancy without altering the original timing measurements.
