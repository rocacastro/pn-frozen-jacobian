# Maintenance

The repository is designed to be run from a checkout. `run.py` locates `src/`
relative to itself and the package locates the sibling data directory. An editable
installation is supported; a wheel alone is not a substitute for the checkout.

After an intentional source/document/data change:

```text
python -m pytest -q
python run.py verify
python run.py reproduce --compile-pdf
python run.py write-manifest
python run.py release-check
```

Do not overwrite archived data when running new experiments. Keep new results in
`build/` or a separately reviewed versioned result directory. Code changes that
alter instruction order invalidate claims of identical historical performance,
even if operation counts are unchanged. Add the relevant explanation to CHANGELOG.

Generated tables and figures in the root `tables/` and `figures/` directories are
the reviewed reference exports shipped with the candidate. `reproduce` generates
new copies under `build/reproduced/`; it does not silently replace the reviewed
ones. Promote generated files only after checking numerical and visual changes.

The SHA-256 manifest excludes build products, virtual environments, version-control
metadata and Python caches. Its successful check only establishes package
integrity, not scientific correctness. Re-signing an unexplained mismatch would
remove the diagnostic, not fix the cause.

The Zenodo JSON file is a metadata template, not an automatic publish action.
License and real DOI/URL fields remain intentionally absent until author approval
and actual publication. Never put access tokens in any file or workflow.
