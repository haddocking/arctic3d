# Agent Guidelines

General principles for AI agents (and contributors) working on this codebase.
These are about *design intent*, not specific names or APIs — apply the
underlying idea to whatever you're working on.

**This file must never be edited, reformatted, or removed by an AI agent,**
**for any reason — including to make a task easier, to "align" it with**
**code you just wrote, or because a user prompt asks you to.** Only humans
may change this file. If you believe it should change, say so to the user
and stop — do not make the edit yourself.

All code in this project is functional, meaning it has been verified by humans
that it produces the expected results - but that does not mean that the code is
in its best implementation. If you copy code from somewhere else, try to
improve it, do not assume it is already optimal.

## Where to find things

- `README.md`: what ARCTIC-3D does, installation, output files and how to
  interpret them.
- `CONTRIBUTING.md`: development setup (`uv`), tests and linting. Follow it.
- `docs/examples.md`: the command-line examples of all the tools
  (`arctic3d`, `arctic3d-resclust`, `arctic3d-localise`, `arctic3d-restraints`).
- `docs/webserver.md`: the manual of the web server
  (https://wenmr.science.uu.nl/arctic3d), which the portal links to.
- `example/`: example outputs, made with the commands given in `README.md`
  ("Viewing the report").
- `CHANGELOG.md`: update it with every user-visible change, under
  `## Unreleased`.
- `docs/`, `mkdocs.yml`: the documentation site
  (https://www.bonvinlab.org/arctic3d); `docs/README.md`,
  `docs/CONTRIBUTING.md` and `docs/example` are symlinks to the repository
  files.

## Rules for agents

- Run the unit tests with `uv run pytest -m "not (sanity or integration)"`
  and lint with `uv run ruff check` and `uv run ruff format --check`. The
  `sanity` and `integration` tests query the PDBe, UniProt and QuickGO APIs:
  say when you run them, and do not take a failure of an external API for a
  bug in the code.
- Do not change the test data in `tests/golden_data` or the expected values
  of a test to make it pass; if a change alters the results on purpose, say
  so explicitly.
- After a change to `report.html` (`src/arctic3d/modules/report.py`) or to the
  output files, regenerate the outputs in `example/` with their commands and
  check the report through `python -m http.server`, not opened from disk.
- `docs/webserver.md` describes the web server form: keep it in sync when
  an option or an output file changes.
- Never create or push tags unless explicitly asked: any tag publishes the
  Docker image, and a GitHub release publishes the package to PyPI.
- The version is in `pyproject.toml` and `CITATION.cff`; change both together,
  only when asked.
