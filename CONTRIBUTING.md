# Developing ARCTIC-3D

## Installation

We use `poetry` to manage the dependencies and the virtual environment, so it makes things easier if you need to install it first; check the [official documentation](https://python-poetry.org/docs/#installation) for more details.

Clone the repository and install the dependencies:

```text
git clone https://github.com/haddocking/arctic3d.git && cd arctic3d
poetry install
```

OR if you prefer `pip`

```text
git clone https://github.com/haddocking/arctic3d.git && cd arctic3d
pip install . pytest-cov
```

## Testing

```text
pytest --cov=./ --cov-report=xml -v
```

## Linting

We use `trunk` as the "all-purpose" linting tool, check its [documentation](https://docs.trunk.io/docs/install).

To check for code style issues, run:

```text
trunk check
```

To automatically fix the issues, run:

```text
trunk fmt
```

If you are using VSCode, then [this extension](https://marketplace.visualstudio.com/items?itemName=Trunk.io) make it easy to check for style errors.
