# Maintainer release checklist

For a future versioned release:

- [ ] Run `release-check`, `verify`, the automated tests, and the relevant data validators.
- [ ] Confirm GitHub Actions pass on all configured platforms.
- [ ] Review documented known issues and preserve archived measurements unchanged.
- [ ] Confirm reproducibility limits for unavailable historical drivers or raw samples.
- [ ] Inspect the computational supplement and figure labels.
- [ ] Review `CITATION.cff`, licensing metadata, and release notes.
- [ ] Create a new semantic-version tag; never move or overwrite an existing tag.
- [ ] Confirm that Zenodo archives the new release as a new version.
- [ ] Record the new version-specific DOI only after Zenodo has minted it.
- [ ] Keep all existing tags and version-specific DOIs immutable, including `v1.0.0` / `10.5281/zenodo.22950468` and `v1.1.0` / `10.5281/zenodo.23090740`.
