"""Version information."""

import re
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


def _read_version() -> str:
    try:
        return version("arctic3d")
    except PackageNotFoundError:
        pyproject_path = Path(__file__).resolve().parents[2] / "pyproject.toml"
        pyproject_contents = pyproject_path.read_text(encoding="utf-8")
        match = re.search(r'^version = "([^"]+)"$', pyproject_contents, re.MULTILINE)
        if match is None:
            raise RuntimeError(f"Could not determine version from {pyproject_path}.")
        return match.group(1)


VERSION = _read_version()
version_match = re.match(r"^(\d+)\.(\d+)\.(\d+)(?:[-+].*)?$", VERSION)
if version_match is None:
    raise RuntimeError(f"Could not parse semantic version from {VERSION!r}.")

v_major, v_minor, v_patch = version_match.groups()
