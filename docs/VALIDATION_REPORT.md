# Validation report — manuscript v0.41 / software 1.0.0-rc2

## Executed checks

| Check | Observed result | Scope |
|---|---:|---|
| Archived data translation | 46 CSV datasets; 23,302 numeric strings preserved | Verified against source-derived digests |
| Scalar model verification | 1,682 checks passed | Costs, reference orders, optima and thresholds |
| Canonical matrices | 15/15 fingerprints match | Five fields, three matrices each |
| Full practical trajectories | 606/606 configurations match | Cycles, convergence and available resource counters |
| Stationary macrocycle counters | 16/16 match | Historical inherited states are not claimed recovered |
| High-precision COC | 16/16 rows within archived rounding | Both English programs actually executed |
| Dense structural audit | 15 configurations recomputed | Full ranks and analytic-Jacobian checks |
| Complex-step relative error | maximum 1.309e-16 (rounded upward) | This validation machine |
| Automated tests | 67 passed | Model identities, domains, solvers, archives, exports and raw-sample smoke test |
| Fresh runner smoke | Eight suites executed | 55 method configurations, 165 samples; not article timing evidence |
| English reproduction | 22 CSV table exports; 13 PDF/PNG panels | Archived values plus analytical model |
| English supplement | Three pages, compiled and visually checked | Three tables and two logical figures |

The corresponding machine-readable evidence is in `metadata/validation/`.
The local system was Linux, Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0,
Matplotlib 3.10.8, mpmath 1.3.0 and pytest 9.0.2. Trajectory, structural and
smoke tests used an explicit one-thread BLAS limit. These are validation settings,
not claims about the user's historical Windows timings. The environment record
also reports the library configuration outside such temporary limits.

## Deliberately non-passing archived-work comparison

A separate audit checked 2,344 archived work values against the current formula.
320 values in the two G5 fusion result files differ, for the reason documented in
`KNOWN_ISSUES.md`. This is not included among the 1,682 passing model checks.
The report is emitted even though the reference data are intentionally preserved.

No elapsed-time sample has been invented or replaced. No new run here is used as
an update to a historical table. The 67 tests include an explicit test that this
known discrepancy remains visible and limited to the identified datasets.

## What these tests do not establish

They do not prove the convergence theorems, historical source-code identity,
bitwise agreement of all iterates, recovery of missing raw repetition timings,
or a universal performance ranking. The historical third-party vector-audit
ledger is not fully regenerated. GitHub Actions is configured but has not been
executed remotely for a repository that has not yet been uploaded.

Author-side Windows verification, review of the work convention, licensing and
actual public deposition are still required before final archival release.
