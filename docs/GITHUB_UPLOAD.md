# Uploading the reviewed candidate to GitHub

The repository has not been published automatically. This guide does not presume
a GitHub username, repository URL, access token, or DOI.

1. Extract the ZIP and enter its `pn-frozen-jacobian-v041` directory.
2. Run the quick-start integrity, verification, and test commands in README.md.
3. Read KNOWN_ISSUES.md and RELEASE_CHECKLIST.md. Select the licenses before the
   public release. Keep original data and any approved corrections distinguishable.
4. Create an empty GitHub repository under the intended owner. Do not add a remote
   README that would conflict with the local one. Copy its actual HTTPS URL.

With Git installed, run the following in the extracted directory. Replace the
placeholder with the URL actually returned by GitHub; do not leave it literal.

```text
git init
git add .
git commit -m "Add English computational companion for manuscript v0.41"
git branch -M main
git remote add origin REPLACE_WITH_ACTUAL_GITHUB_REPOSITORY_URL
git push -u origin main
```

No `.git` folder or credentials are shipped. The included workflow tests the
repository on Linux and Windows after upload; the local build cannot claim that
remote GitHub Actions jobs have already passed.

After any intentional edits to tracked files, re-run the tests and update the
manifest with `python run.py write-manifest` before the next commit. Never use
that command to dismiss an unexplained integrity mismatch.

Once the author has validated the candidate on the experimental computer and
approved the work-convention issue, create a real release/tag. Connect or deposit
the repository in Zenodo using the actual release. Only then insert the real URL,
DOI, date and final version in CITATION.cff and the manuscript reference. No
placeholder identifier should be represented as an existing publication.
