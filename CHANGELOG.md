# Changelog

## Unreleased

- Publish the documentation at bonvinlab.org/arctic3d with MkDocs - PR #534

## 0.7.0

- New report.html with a Mol* view of the binding surfaces and the results
  (`--no_report` to skip); new dependency: molviewspec - PR #523
- `arctic3d.log` has a header with the version and a footer with the run
  time - PR #523
- Python 3.9 is no longer supported - PR #523
- report.html shows the `--biological_clustering` results - PR #526
- Fixed `--biological_clustering` with several types - PR #525
- Fixed black bars in sequence_probability.html - PR #527
- Fixed failing bar plot images: matplotlib replaces kaleido - PR #524
- Example output with `--biological_clustering` - PR #529
- Documentation updated; code linted with ruff - PR #531

## 0.6.0

- New `--biological_clustering` option to cluster partners by one or more
  biological attributes: `location` (subcellular location), `function`
  (protein function) and/or `process` (biological process) - PR #510
- Guidance in the docs on how to read and interpret the ARCTIC-3D output -
  PR #511
- The CLI help shows the default parameter values - PR #509
- `arctic3d_resclust` can write a plain-text summary of the clustering
  results - PR #343
- Fixed the remote BLAST query, and updated for the deprecated Plotly
  axis-title API - PR #512
- Corrected chain handling - PR #415
- Added a dedicated `db` directory and improved help messages - PR #402
- Updated BLAST installation instructions and handling - PR #394
- Removed Poetry in favor of a simpler build setup - PR #493
- GitHub Actions workflows to publish to Test PyPI and PyPI - PR #407
- Moved URLs into a dedicated module (PR #428) and refactored functions that
  relied heavily on `None` (PR #429)
- Removed `defusedxml` (PR #424) and redundant download-skipping logic
  (PR #423)
- `stale.yml` workflow for stale issues and PRs, and general CI updates -
  PR #504
- Dependency updates, including Plotly 5 → 6, Biopython 1.83 → 1.85,
  pandas 2.2.3, jsonpickle 4.x, lxml 5.4, and newer matplotlib, scipy,
  requests and test tooling

## 0.5.1 and earlier

See the [GitHub releases](https://github.com/haddocking/arctic3d/releases).
