import re
from pathlib import Path

from arctic3d.version import VERSION, v_major, v_minor, v_patch


def test_version_matches_pyproject():
    pyproject_path = Path(__file__).resolve().parents[1] / "pyproject.toml"
    pyproject_contents = pyproject_path.read_text(encoding="utf-8")
    expected_version = re.search(
        r'^version = "([^"]+)"$', pyproject_contents, re.MULTILINE
    ).group(1)

    assert VERSION == expected_version
    assert (v_major, v_minor, v_patch) == tuple(expected_version.split("."))
