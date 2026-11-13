# Developing ARCTIC-3D

## Installation

We use `uv` to manage the dependencies and the virtual environment, so it makes things easier if you need to install it first; check the [official documentation](https://docs.astral.sh/uv/getting-started/installation/) for more details.

Clone the repository and install the dependencies:

```text
git clone https://github.com/haddocking/arctic3d.git && cd arctic3d
uv sync
```

OR if you prefer `pip`

```text
git clone https://github.com/haddocking/arctic3d.git && cd arctic3d
pip install .
pip install pytest pytest-cov pytest-mock coverage hypothesis mypy ruff types-requests
```

## Testing

```text
uv run pytest -m "not (sanity or integration)"
```

The `sanity` (API checks) and `integration` tests need an internet connection: run them with `uv run pytest -m sanity` and `uv run pytest -m integration`.

## Linting

The code is linted and formatted with [ruff](https://docs.astral.sh/ruff/), configured in `pyproject.toml` (lines of at most 79 characters). The CI runs both checks:

```text
uv run ruff check        # add --fix to fix what can be fixed automatically
uv run ruff format       # or ruff format --check to only check
```
