from pathlib import Path

import pytest

from arctic3d.cli_resclust import main

from . import golden_data


@pytest.fixture
def example_pdbpath():
    """Example pdb path."""
    return Path(golden_data, "1rypB_r_b.pdb")


def test_wrong_calls(example_pdbpath):
    main(example_pdbpath, "100,101,102,133,134,135", None, None)
