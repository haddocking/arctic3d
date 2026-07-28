"""Version information."""

import re
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


def _find_pyproject_path() -> Path:
    for directory in Path(__file__).resolve().parents:
        pyproject_path = directory / "pyproject.toml"
        if pyproject_path.exists():
            return pyproject_path

    raise RuntimeError("Could not find pyproject.toml for version fallback.")


def _read_version() -> str:
    try:
        return version("arctic3d")
    except PackageNotFoundError:
        pyproject_path = _find_pyproject_path()
        pyproject_contents = pyproject_path.read_text(encoding="utf-8")
        match = re.search(r'^version = "([^"]+)"$', pyproject_contents, re.MULTILINE)
        if match is None:
            raise RuntimeError(f"Version field not found or malformed in {pyproject_path}.")
        return match.group(1)


VERSION = _read_version()
_semver_match = re.match(r"^(\d+)\.(\d+)\.(\d+)(?:[-+].*)?$", VERSION)
if _semver_match is None:
    raise RuntimeError(f"Could not parse semantic version from {VERSION!r}.")

v_major: str
v_minor: str
v_patch: str
v_major = _semver_match.group(1)
v_minor = _semver_match.group(2)
v_patch = _semver_match.group(3)
