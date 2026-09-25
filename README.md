# Predictor-memory frozen-Jacobian methods

**Computational repository for _Frozen-Jacobian methods with predictor and memory: order, optimal selection, and fusion_**

Rodrigo Castro Marín · Instituto de Matemáticas, Universidad de Valparaíso, Chile

This repository contains implementations, archived numerical data, reproduction
scripts, figures, and supplementary computational material supporting the study
of predictor-memory frozen-Jacobian methods for nonlinear systems.

The principal methods are $P_N$, Shamanskii's family $S_m$, the fused family
$\widehat P_{N,q}$, and the external comparators $M6$ and $M8$. Throughout the
repository, $q$ denotes the number of terms in the truncated Neumann sum, so the
corresponding polynomial degree is $q-1$.

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
Activation of the virtual environment is optional. These commands may be run
directly from a checkout; alternatively, use an editable installation with
`python -m pip install -e .`.

`reproduce` reads the archived data and regenerates repository tables, figures,
and supplementary LaTeX. It does **not** repeat the publication timing
experiments. Generated files are written to `build/reproduced/`. The compiled
supplement is available in [`supplement/`](supplement/). If `pdflatex` is
installed, rebuild its PDF with:

```powershell
.venv\Scripts\python run.py reproduce --compile-pdf
```

## Reproduction, verification, and new measurements

These are deliberately separate tasks.

| Task | Command | Purpose |
|---|---|---|
| Verify archived values and algorithms | `python run.py verify` | Check numerical integrity, cost formulas, canonical matrices, stopping rules, and resource counts |
| Rebuild tables and figures | `python run.py reproduce` | Regenerate published tables, figures, and supplementary material from archived data |
| Make new timing measurements | `python run.py benchmark --suite all` | Produce new timing samples on the current computer and store them separately |

For a quick execution test, not a publication benchmark:

```text
python run.py benchmark --suite all --quick --threads 1 --outdir build/smoke
```

For a full fresh timing run:

```text
python run.py benchmark --suite all --outdir build/current
```

The default does not alter the BLAS thread setting. Use `--threads N` only as an
explicit experimental choice and record it with the resulting environment data.
The archived Windows environment files report 16 configured OpenBLAS threads.
New timing measurements must be kept separate from the archived tables.

Additional commands are available for deterministic analyses:

```text
python run.py models
python run.py structure
python run.py orders
```

`models` recomputes parameter selection and affine interval checks. `structure`
recomputes the dense-field structural audit. `orders` runs the separate
high-precision order-verification protocols and checks the resulting COC values
against the archived rounded values.

## Repository validation

The repository validation suite checks the following archived material:

- 46 numerical CSV datasets, with 23,302 original numeric cell strings preserved;
- 606 archived configurations reproduced in cycle counts, convergence status,
  and available resource counters;
- 15 canonical matrix fingerprints and 1,682 scalar cost-model checks;
- 16 stationary resource configurations checked independently of timing data;
- 16 high-precision COC rows reproduced within the archived rounding accuracy.

Machine-readable validation results are stored in
[`metadata/validation/`](metadata/validation/). The scope and limitations of
these checks are described in
[`docs/VALIDATION_REPORT.md`](docs/VALIDATION_REPORT.md).

## Provenance and known issue

Read [`docs/KNOWN_ISSUES.md`](docs/KNOWN_ISSUES.md) before using the archived
fusion work columns. In two archived fusion datasets, the $G_5$ Jacobian
evaluation cost was recorded as

\[
2+\frac{18.2}{n},
\]

whereas the canonical cost model used in the article is

\[
2+\frac{17+\kappa}{n}.
\]

The archived values are preserved for provenance, and the repository provides
separately recomputed canonical costs. No elapsed-time measurement is altered by
this discrepancy.

The archived data do not contain every original historical timing driver or the
individual repetition samples for every experiment. The float64 implementation
provided here is therefore a documented reimplementation checked against the
available deterministic outputs. The initialization summary is retained from
archived summary values because the individual timing samples for that block are
not available. The full reproducibility scope is described in
[`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md).

## Repository layout

```text
src/pn_fusion/       Methods, cost models, validation, and report generators
experiments/         High-precision order-verification programs
run.py               Command-line entry point
configs/             Fixed experimental settings
tests/              Automated tests
data/reference/     Archived numerical data; never overwritten by experiments
data/manuscript/    Values transcribed from article tables when raw samples are unavailable
tables/             CSV exports for article tables
figures/            PDF and PNG figures
supplement/         Supplementary PDF, LaTeX source, and figure files
metadata/           Provenance, validation records, environments, and SHA-256 manifest
docs/               Methods, protocols, data dictionary, manuscript map, and release notes
```

Public repository file names, code comments, command-line interfaces, plot
labels, and documentation are in English. Mathematical identifiers such as
`H3`, `Phat2`, and `mu1` retain their definitions from the article.

## Supplementary material

The reader-facing supplementary document is available as
[`supplement/supplementary_computational_material.pdf`](supplement/supplementary_computational_material.pdf).
It contains detailed timing tables, the full sensitivity table for the optimal
family members, and secondary graphical comparisons omitted from the main
article for concision. The underlying full-precision data remain available in
`data/reference/`.

## Citation

Software citation metadata are provided in [`CITATION.cff`](CITATION.cff).
After archival of the public release in Zenodo, the DOI will be added to the
repository metadata.

For licensing information, see [`LICENSE_POLICY.md`](LICENSE_POLICY.md).
