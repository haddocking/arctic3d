from pathlib import Path

import pytest

import os
import shutil

from arctic3d.cli_resclust import main

from . import golden_data


@pytest.fixture
def example_pdbpath():
    """Example pdb path."""
    return Path(golden_data, "1rypB_r_b.pdb")


def test_resclust_cli(example_pdbpath):
    main(
        example_pdbpath,
        "100,101,102,133,134,135",
        None,
        7.0,
        "average",
        "distance",
        None,
    )


def test_wrong_residue_list(example_pdbpath):
    with pytest.raises(SystemExit) as e:
        main(
            example_pdbpath,
            "100,101,102,133,134,135%",
            None,
            9.0,
            "average",
            "distance",
            None,
        )
    assert e.type == SystemExit
    assert e.value.code == 1


def test_resclust_maxclust(example_pdbpath):
    main(
        example_pdbpath,
        "100,101,102,133,134,135",
        None,
        2,
        "average",
        "maxclust",
        None,
    )


def test_resclust_genoutput(example_pdbpath):
    main(
        example_pdbpath,
        "100,101,102,133,134,135",
        None,
        2,
        "average",
        "maxclust",
        "resclustout",
    )
    assert os.path.exists("resclustout") == True
    assert os.path.exists("resclustout/clustered_residues.out") == True
    shutil.rmtree("resclustout")

