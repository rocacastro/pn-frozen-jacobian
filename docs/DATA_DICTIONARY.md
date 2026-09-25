# Data dictionary and units

`data/reference/` contains archived scientific tables, not a claim of individual
raw experimental samples. All files use UTF-8 CSV with a header. The full
per-file header list is `metadata/dataset_schema.json`. Provenance and numerical
integrity are recorded in `metadata/data_provenance.json`.

| Column or convention | Meaning |
|---|---|
| `system`, or `family` in fusion tables | Field identifier H1--H5, G1/G5, or Hammerstein |
| `family` in the family sweep | P or S, not a vector-field name |
| `N`, `index`, `m` | Fixed method correction count; interpret together with method/family |
| `q` | Number of Neumann terms; polynomial degree q-1 |
| `cycles` | Ordinary outer cycles, except fused/composed methods where it counts macrocycles |
| `macrocycles` | Explicit macrocycle count in the archived fusion tables |
| `F_total_including_terminal` | All distinct F evaluations, including the last stopping test |
| `F` in stationary resource data | F evaluations inside a stationary macrocycle; no terminal stopping test |
| `J`, `LU`, `solves`, `matvecs` | Complete evaluations/factorizations/solves/matrix-vector products |
| `matvec_entries` | Number of charged scalar products, n squared per dense matvec |
| `matrix_scale_entries`, `vector_scale_entries` | Charged scalar entries, not matrix/vector operation calls |
| `residual`, `error` | Euclidean norms, respectively of F(x) and x-alpha |
| `time_median_s`, `time_iqr_s`, `time_min_s` | Elapsed time in seconds |
| `*_ms` | Milliseconds; the manuscript tables often use these instead of seconds |
| `work*`, `cost*` | Equivalent products under the indicated kappa, not seconds |
| `eta*` | Natural-log reference index log(p_ref)/cost |
| `COC*` | Measured log-error-ratio convergence estimate, distinct from p_ref |
| `clear*`, `*_2iqr` | Descriptive separation flag under the manuscript's two-IQR rule |

A suffix `_pct` explicitly denotes a percentage (100 times a ratio minus one).
Some historical columns such as `CPU_relative_P5_vs_M6`,
`work_relative_P5_vs_M6`, and `relative_eta_advantage_P` are **fractions**, not
percentages. Their plot/table exporters multiply by 100. Always inspect the
export function rather than assuming all columns named `relative` use one unit.
Pair orientation A/B is preserved: negative time/work differences favor A;
positive index differences favor A.

Missing cells, `--`, and NaN-style diagnostics retain their source meaning; they
are not newly imputed measurements. Textual method labels may retain compact
mathematical spellings such as `hatP2` or `P2oP2`. Decimal strings have not been
rounded during translation. Full filenames and checksums identify the source
unambiguously even where a column name is shared by different experiments.

The G5 fusion work discrepancy is documented separately. Reproduction exports
preserve the archived work columns; canonical recomputation is not mixed into
the same column without a distinct name and provenance.
