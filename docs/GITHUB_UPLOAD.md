# Maintainer GitHub and Zenodo guide

The public repository is `https://github.com/rocacastro/pn-frozen-jacobian`.

The initial `v1.0.0` release is archived in Zenodo under DOI
`10.5281/zenodo.22950468`.

For a subsequent public release:

1. Run the integrity, verification, test, and reproduction commands documented in the repository.
2. Review `KNOWN_ISSUES.md` and `RELEASE_CHECKLIST.md`.
3. Confirm GitHub Actions pass on all configured platforms.
4. Review the release diff and ensure new data are accompanied by their protocols, environment records, and validation metadata.
5. Create a new semantic-version tag without modifying any existing tag.
6. Publish the GitHub release and confirm that Zenodo creates a new version of the existing record.
7. Add the newly minted version-specific DOI to release metadata only after it is verified.

After an intentional edit to tracked files, regenerate the integrity manifest with
`python run.py write-manifest` only after reviewing the change. Never regenerate
the manifest merely to suppress an unexplained mismatch.
