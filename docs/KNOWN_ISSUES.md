# Known issues and non-silent corrections

## K1. Archived G5 Jacobian evaluation cost

**Status: disclosed; archived numerical values are preserved unchanged.**

The article and the canonical implementation specify

$$\mu_1(G_5;\kappa)=2+\frac{17+\kappa}{n}.$$

The archived files `fusion_stationary_results.csv` and
`fusion_end_to_end_results.csv` use $\mu_1=2+18.2/n$ independently of kappa in
their work columns. Direct reconstruction from the stored F, J, LU, solve and
matvec counts identifies this difference; no inference from noisy timings is
needed. For a stored row with two J evaluations,

$$W_{archived}-W_{canonical}=2n(1.2-\kappa).$$

At kappa=1 this is +0.4n: +16, +22, +40 and +80 equivalent products for
n=40,55,100 and200 respectively. There are 320 discrepant work values in these
two files across the archived kappa controls, out of 2,344 checked work values.
The maximum absolute relative total-work discrepancy is 0.184086% across all
those controls; at kappa=1 it is 0.018993%.

Both methods in each affected comparison evaluate two Jacobians. Thus the
**absolute cost difference and its sign are unchanged**, as is the derived
stationary threshold. Their elapsed times, iteration counts and operation
counts are unchanged. Relative-work percentages can change because their
denominator changes. At kappa=1 the largest change is
0.001686 percentage points.
All 40 affected pair signs and their percentages rounded to two decimal places
agree at kappa=1. Higher-precision values differ; do not claim complete numerical
identity of the relative work tables across all conventions.

### Files and behavior

- `data/reference/` retains every archived work value.
- `models.py` and fresh benchmark work use the canonical formula.
- `python run.py verify` writes a separate work discrepancy ledger instead of
  rewriting the archive or pretending that all work cells pass.
- `metadata/validation/archived_work_discrepancies.csv` lists every differing value.
- `metadata/validation/fusion_work_comparison_kappa1.csv` shows old and canonical
  relative percentages and tests their signs and two-decimal rounding.
- Figures/tables labeled as archived reproduction continue to use archived work.

Any corrected derived table should be distributed as a separate, clearly
identified artifact with an explicit provenance link; the archived source data
remain unchanged.

## K2. Historical code and repetition samples

Some original performance drivers and individual repetition timings were not
present in the supplied package. Deterministic equivalence of the new English
implementation has been checked, but its timer overhead and historical instruction
identity are not established. The original stationary inherited matrices are also
absent. See REPRODUCIBILITY.md for exactly what can and cannot be regenerated.

## K3. Archival scope

The initial public software release is archived under DOI
`10.5281/zenodo.22950468`. Data added after that release are not part of that
version-specific archive and are identified separately until a subsequent
versioned release is created.
