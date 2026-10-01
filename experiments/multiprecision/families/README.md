# Multiprecision family experiment

This directory contains the high-precision implementation used for the predictor-memory/Shamanskii family comparison.

The distributed dataset uses:

- `H1` (`n=20`): `P3` versus `S3`;
- `H5` (`n=200`): `P9` versus `S10`;
- delayed initialization `D`;
- starting distance `0.30`;
- tolerances `1e-50`, `1e-200`, `1e-800`, `1e-1024`;
- 1300 working digits and 2600 verification digits;
- two warmups and 31 paired repetitions per group.

Install and run the self-test with:

```console
python -m pip install -r requirements.txt
python run_family_multiprecision.py test
```

A fresh experiment should use a new output directory. The archived measurements are under `data/multiprecision/families/` and must not be overwritten by new runs.
