# Maintainer GitHub and Zenodo guide

The public repository is `https://github.com/rocacastro/pn-frozen-jacobian`.

Before an archival release:

1. Run the integrity, verification, test, and reproduction commands in `README.md`.
2. Review `KNOWN_ISSUES.md` and `RELEASE_CHECKLIST.md`.
3. Select and add the final licenses.
4. Confirm GitHub Actions pass on all configured platforms.
5. Enable this repository in the Zenodo GitHub integration.
6. Create the reviewed GitHub tag and release `v1.0.0`.
7. Confirm that Zenodo archives the release and assigns a DOI.
8. Add the verified DOI to `CITATION.cff`, the README, and the article citation.

After any intentional edit to tracked files, regenerate the integrity manifest with
`python run.py write-manifest` only after reviewing the change. Never regenerate the
manifest merely to suppress an unexplained mismatch.
