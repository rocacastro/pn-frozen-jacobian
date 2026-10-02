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

The supplement archived with the initial public release remains available as
[`supplement/supplementary_computational_material.pdf`](supplement/supplementary_computational_material.pdf).
The current computational supplement, synchronized with the submission-ready article, is available at
[`supplement/current/Supplementary_Computational_Material.pdf`](supplement/current/Supplementary_Computational_Material.pdf).
Underlying archived and multiprecision data are retained in their documented
data directories.

## Citation

Software citation metadata are provided in [`CITATION.cff`](CITATION.cff).
The initial public release is archived in Zenodo under DOI `10.5281/zenodo.22950468`. The multiprecision datasets were first archived in version `v1.1.0` under DOI `10.5281/zenodo.23090740`. The archive series is available through the Zenodo concept DOI `10.5281/zenodo.22950467`.

For licensing information, see [`LICENSE_POLICY.md`](LICENSE_POLICY.md).

## Multiprecision datasets

The repository includes complete high-precision timing datasets for the independently selected predictor-memory/Shamanskii family members on `H1` and `H5`, and for the external pairs `P5/M6` and `P7/M8` on the same fields. The archived data include individual paired timings, higher-precision controls, fixed protocols, session/platform metadata, and descriptive audit summaries.

Start with [`docs/MULTIPRECISION_DATA.md`](docs/MULTIPRECISION_DATA.md) and [`docs/REPRODUCIBILITY_MULTIPRECISION.md`](docs/REPRODUCIBILITY_MULTIPRECISION.md). Validate the added files with:

```console
python tools/validate_multiprecision_data.py
```

The multiprecision datasets were first archived in `v1.1.0` under DOI `10.5281/zenodo.23090740` and are preserved unchanged in `v1.1.1`. The `v1.1.1` patch synchronizes the current computational supplement and public documentation; it does not change the multiprecision measurements or protocols.
