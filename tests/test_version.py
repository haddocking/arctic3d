import re
from pathlib import Path

from arctic3d.version import VERSION, v_major, v_minor, v_patch


def test_version_matches_pyproject():
    pyproject_path = Path(__file__).resolve().parents[1] / "pyproject.toml"
    pyproject_contents = pyproject_path.read_text(encoding="utf-8")
    version_match = re.search(
        r'^version = "([^"]+)"$', pyproject_contents, re.MULTILINE
    )
    assert version_match is not None

    expected_version = version_match.group(1)
    semver_match = re.match(r"^(\d+)\.(\d+)\.(\d+)(?:[-+].*)?$", expected_version)
    assert semver_match is not None

    assert VERSION == expected_version
    assert all(isinstance(part, str) for part in (v_major, v_minor, v_patch))
    assert (v_major, v_minor, v_patch) == semver_match.groups()
