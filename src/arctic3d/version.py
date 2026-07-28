"""Version information."""

import re
from functools import lru_cache
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


@lru_cache(maxsize=1)
def _find_pyproject_path() -> Path:
    start_path = Path(__file__).resolve()
    for directory in start_path.parents:
        pyproject_path = directory / "pyproject.toml"
        if pyproject_path.exists():
            return pyproject_path

    raise RuntimeError(f"Could not find pyproject.toml for version fallback from {start_path}.")


def _parse_version_parts(version_string: str) -> tuple[str, str, str]:
    version_parts_match = re.match(
        r"^(\d+)\.(\d+)\.(\d+)(?:[-+].*)?$", version_string
    )
    if version_parts_match is None:
        raise RuntimeError(
            "Version field exists but does not match expected format "
            f"'X.Y.Z[...suffix]': {version_string!r}."
        )

    return (
        version_parts_match.group(1),
        version_parts_match.group(2),
        version_parts_match.group(3),
    )


def _read_version() -> str:
    try:
        return version("arctic3d")
    except PackageNotFoundError:
        pyproject_path = _find_pyproject_path()
        pyproject_contents = pyproject_path.read_text(encoding="utf-8")
        match = re.search(r'^version = "([^"]+)"$', pyproject_contents, re.MULTILINE)
        if match is None:
            if 'version = "' in pyproject_contents:
                raise RuntimeError(f"Version field is malformed in {pyproject_path}.")
            raise RuntimeError(f"Version field not found in {pyproject_path}.")
        return match.group(1)


VERSION = _read_version()
v_major: str
v_minor: str
v_patch: str
v_major, v_minor, v_patch = _parse_version_parts(VERSION)
