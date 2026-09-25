# Predictor-memory frozen-Jacobian methods

**Computational companion to manuscript v0.41 · local release candidate 1.0.0-rc2**

Rodrigo Castro Marín · Instituto de Matemáticas, Universidad de Valparaíso, Chile

This repository provides English implementations, archived numerical data,
reproduction scripts, figures, and supplementary computational material for
*Frozen-Jacobian methods with predictor and memory: order, optimal selection, and
fusion* (English rendering of the title of manuscript v0.41).

The core methods are $P_N$, Shamanskii $S_m$, the fused family
$\widehat P_{N,q}$, and the external comparators $M6$ and $M8$.
Throughout, $q$ is the number of Neumann terms; the polynomial degree is $q-1$.
The manuscript is maintained separately. This package does **not** silently
translate or replace the current Spanish article.

## Start here

Use Python 3.11 or later. From the repository root, on Windows:

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
.venv\Scripts\python run.py release-check
.venv\Scripts\python run.py verify
.venv\Scripts\python -m pytest -q
.venv\Scripts\python run.py reproduce
```

On Linux or macOS, replace `.venv\Scripts\python` with `.venv/bin/python`.
Activation of the virtual environment is optional; these commands do not require
changing PowerShell's execution policy. Run from the checkout, or use an editable
installation (`python -m pip install -e .`); this is a research checkout with data,
not a standalone wheel containing the archived results.

`reproduce` reads the archived data; it does **not** remeasure publication timings.
It writes English CSV table exports, PDF/PNG figure panels, and supplementary
LaTeX to `build/reproduced/`. The already compiled English supplement is in
[`supplement/`](supplement/). With `pdflatex` installed, rebuild its PDF with:

```powershell
.venv\Scripts\python run.py reproduce --compile-pdf
```

## Three different tasks

| Task | Command | What it establishes |
|---|---|---|
| Verify archived values and algorithms | `python run.py verify` | Numeric integrity, cost formulas, canonical matrices, stopping and resource counts |
| Rebuild the published evidence | `python run.py reproduce` | Tables and English figures derived from archived aggregate data |
| Make new timing measurements | `python run.py benchmark --suite all` | New individual samples on the current computer, saved separately |

For a quick execution test, not a publication benchmark:

```text
python run.py benchmark --suite all --quick --threads 1 --outdir build/smoke
```

For the full fresh timing protocol:

```text
python run.py benchmark --suite all --outdir build/current
```

The default does not change the BLAS thread setting. Specify `--threads N` only
as an explicit experimental choice and record it. The original supplied Windows
environment files report 16 configured OpenBLAS threads. Fresh times from this
reimplementation must not be spliced into the historical tables.

Additional commands:

```text
python run.py models
python run.py structure
python run.py orders
```

`models` recomputes parameter selection and affine interval checks. `structure`
recomputes the dense-field audit. `orders` runs the separate, potentially slower
high-precision protocols and checks their COC values against the archived rounding.

## Validation performed for this candidate

- 46 numerical CSV datasets, with 23,302 original numeric cell strings preserved.
- 606 archived configurations reproduced in cycles, convergence status, and
  available resource counters; this is not a claim of bitwise residual agreement.
- 15 canonical matrix fingerprints and 1,682 scalar model checks passed.
- 16 stationary resource configurations checked independently of historical timing.
- Both high-precision programs were executed; all 16 reported COC rows agree
  within the stated archived precision.

Machine-readable results are in [`metadata/validation/`](metadata/validation/).
See [`docs/VALIDATION_REPORT.md`](docs/VALIDATION_REPORT.md) for the test scope.

## Important provenance and known issue

**Read [`docs/KNOWN_ISSUES.md`](docs/KNOWN_ISSUES.md) before using the fusion work
columns.** In two archived fusion datasets, the $G_5$ Jacobian evaluation cost was
held at $2+18.2/n$, while v0.41 declares $2+(17+\kappa)/n$. The original numbers
are retained; the canonical implementation and a separate discrepancy report are
provided. No elapsed-time value has been changed.

Several historical timing drivers and individual repetition samples were absent
from the supplied manuscript package. The float64 code is therefore an explicit
reimplementation checked against the archived deterministic outputs, not a claim
that the exact original timing source has been recovered. The initialization
summary is transcribed from v0.41, not reconstructed from missing raw timings.
Full limits are stated in [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md).

## Repository layout

```text
src/pn_fusion/       English methods, cost model, data checks, and report generators
experiments/        Faithful English high-precision order programs
run.py              Single command-line entry point
configs/            Fixed experimental settings, recorded for inspection
tests/              Automated tests
data/reference/     Archived numbers; never overwritten by experiments
data/manuscript/    Explicitly identified manuscript transcriptions
tables/             English CSV exports for the 22 manuscript tables
figures/            English PDF and PNG panels for the main and supplementary figures
supplement/         English supplementary PDF, LaTeX, and required figure PDFs
metadata/           Provenance, validation records, environments, and SHA-256 manifest
docs/               Methods, protocols, data dictionary, manuscript map, release guide
```

The original Spanish filenames appear only as exact provenance identifiers in
metadata; public file names, code comments, interfaces, plot labels, and documents
are in English. Mathematical labels such as `H3`, `Phat2`, and `mu1` retain their
meaning. Editorial drafts and referee correspondence are not part of this repository.

## Citation and release

[`CITATION.cff`](CITATION.cff) contains the known author and software metadata.
No GitHub URL or Zenodo DOI has been fabricated. The author must select the
license, review the known issue, validate the candidate on the experimental
computer, and then publish the actual repository/deposit identifiers. See
[`docs/GITHUB_UPLOAD.md`](docs/GITHUB_UPLOAD.md),
[`docs/RELEASE_CHECKLIST.md`](docs/RELEASE_CHECKLIST.md), and
[`LICENSE_POLICY.md`](LICENSE_POLICY.md).

The reproducibility objection is **not** represented as closed merely because
this local package has been built.
