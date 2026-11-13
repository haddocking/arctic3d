# Changelog

## Latest changes

- Added report.html generation to visualise results (Mol* view of all
  binding surfaces and of each cluster, residues colored by probability;
  tables of the clusters, residue probabilities and retrieved interfaces; the
  dendrogram and sequence plot; links to the output files); turned off with
  `--no_report`. New dependency: molviewspec
- `arctic3d.log` starts with a header giving the ARCTIC-3D version and ends
  with a footer giving the run time and how to view report.html
- Python 3.9 is no longer supported (molviewspec needs Python >= 3.10)

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
