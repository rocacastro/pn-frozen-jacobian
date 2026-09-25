# Reproducibility scope and provenance

## Canonical source

The repository targets the accompanying 34-page article, with 22 main tables, four main figures, and a three-page computational supplement. `metadata/source_identity.json` records
the hashes of the supplied manuscript and archive. Earlier drafts are not used as
alternative mathematical specifications.

## 1. Reproduction from archived results

`python run.py reproduce` exports the numerical contents of all 22 manuscript
tables and regenerates English plots and supplementary material from archived
values. Table export precision and row layout can differ from the typeset table;
values are not selected to favor a method. Every export has a source manifest.

The 46 CSV files have English filenames and translated textual fields. Numerical
cell strings, including their original exponent formatting, remain unchanged.
`metadata/data_provenance.json` records original filenames, source and translated
file hashes, and numeric-cell digests. This is a translation, not a byte-identical
copy of every original CSV. The exact source names in that manifest are identifiers,
not untranslated documentation.

`data/manuscript/initialization_summary.csv` is a separately labeled transcription
of the twelve initialization rows in the article. Its rounded medians and IQRs are not
raw repetitions. No missing individual samples have been synthesized.

## 2. Deterministic algorithm checks

`python run.py verify` checks 606 archived configurations: cycle counts,
convergence flags, residual acceptance, and available operation counters. It also
checks the 15 canonical matrix fingerprints, 1,682 scalar model values, and 16
stationary resource counts. The result confirms matching aggregate algorithmic
behavior, not equality of every intermediate floating-point bit.

The float64 implementation is reconstructed from the accepted formulas and the
available Hammerstein and Phat2/P10 drivers. Some earlier standalone drivers were
not contained in the supplied archive and were not recovered. The code does not
claim historical identity with those missing programs. Exact instruction order,
Python overhead, timer placement, allocation behavior, and scheduling may differ.

The two high-precision drivers were available. Their English rewrites preserve
their numerical protocols; both were executed for this package. The main order
driver constructs its initial mpmath values before increasing precision, as the
source did. The q-transition driver constructs its decimal starting values after
setting precision. They are intentionally separate datasets, and their distinct
COC estimates must not be merged into a single run.

## 3. Fresh performance experiments

`python run.py benchmark` creates new individual timing samples and summaries in
`build/current/`. It never writes into `data/reference/`. The default protocols
use 31 interleaved repetitions and two warmups; initialization uses 51 repetitions
and three warmups. The fresh shuffle rule is explicit in `docs/PROTOCOLS.md`.
It is not claimed to be the undisclosed historical schedule of missing drivers.

The measurement is elapsed **wall-clock** time from `time.perf_counter_ns`, not
process CPU time. Historical filenames and column names containing `CPU` are
retained as schema identifiers. BLAS threading, backend, hardware, load, thermal
state, and Python/software versions affect the measured values. Matching old
medians or IQRs exactly is neither expected nor used as a test.

For the stationary fusion suite, historical inherited matrices were not archived.
The English runner constructs a valid inherited state with one untimed delayed
macrocycle and then measures a full stationary macrocycle. Resource counts match
the specification. The resulting new times are not reproductions of the old
stationary state or old wall-clock measurements.

## Environment records

The two supplied original environment files describe Hammerstein and Phat2/P10 on
Windows 11, Python 3.13.6, NumPy 2.3.2, SciPy 1.17.1, with OpenBLAS configured for
16 threads. Local absolute paths have been removed. These files are not imputed
to other historical runs without evidence. The repository build/validation
machine has its own separate record in `metadata/validation/environment.json`.

`constraints-reference.txt` pins only NumPy and SciPy versions actually recorded
by those author runs. It is not a full historical environment lockfile. The
validation environment records the additional dependency versions used here.

## Remaining limitations

- Original per-repetition timing samples are absent from the supplied archive.
- The original initialization result CSV was not recovered; the paper summary is
  preserved explicitly as a transcription.
- Several original benchmark runners and stationary inherited states are absent.
- The historical scalar-division microbenchmark and its raw timing samples are
  not supplied here. Its reported kappa values are only documented sensitivity
  controls; the model does not derive wall-clock timings from them.
- `historical_vector_audit.csv` is an archived ledger. The original symbolic
  auditing scripts were absent, so this package does not independently reproduce
  every third-party order claim recorded there.
- The archived G5 work convention is inconsistent with the canonical formula;
  the mismatch is preserved and quantified in `KNOWN_ISSUES.md`.
- The article is maintained separately and is not modified by repository reproduction commands.

These limits are part of the reproducibility record, not failures hidden by
replacing old data with fresh runs. The Zenodo DOI and final license metadata are added only after the archived release is created.
