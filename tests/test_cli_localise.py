import os
import shutil
from pathlib import Path

import pytest

from arctic3d.cli_localise import (
    call_uniprot,
    get_quickgo_information,
    get_uniprot_subcellular_location,
    main,
)

from . import golden_data


@pytest.fixture
def empty_cluster_filepath():
    """Empty cluster filepath."""
    return Path(golden_data, "clustered_interfaces_empty.out")


@pytest.fixture
def example_uniprot_data():
    return {
        "comments": [
            {
                "type": "SUBCELLULAR_LOCATION",
                "locations": [
                    {"location": {"value": "Nucleus", "evidences": []}}
                ],
            }
        ],
        "dbReferences": [
            {
                "type": "GO",
                "id": "GO:0005829",
                "properties": {"term": "C:cytosol", "source": "IDA:HPA"},
            }
        ],
    }


@pytest.mark.integration
def test_localise_cli_empty(empty_cluster_filepath):
    """Test localise cli with empty cluster file."""
    start_cwd = os.getcwd()
    run_dir = "arctic3d-localise"
    main(
        empty_cluster_filepath,
        "arctic3d-localise",
        None,
        None,
        None,
        None,
        None,
    )
    os.chdir(start_cwd)
    # Check that the output directory has been created
    assert Path(run_dir).exists()
    # check existence of log file
    assert Path("arctic3d-localise", "arctic3d-localise.log").exists()
    shutil.rmtree(run_dir)


def test_get_quickgo_information(example_uniprot_data):
    """Test get_quickgo_information."""
    obs_locs = get_quickgo_information(example_uniprot_data, quickgo_key="C")
    exp_locs = ["cytosol"]
    assert exp_locs == obs_locs
    obs_locs_P = get_quickgo_information(example_uniprot_data, quickgo_key="P")
    exp_locs_P = []
    assert exp_locs_P == obs_locs_P


def test_get_uniprot_subcellular_location(example_uniprot_data):
    """Test get_uniprot_subcellular_location"""
    obs_locs = get_uniprot_subcellular_location(example_uniprot_data)
    exp_locs = ["Nucleus"]
    assert exp_locs == obs_locs


def test_call_uniprot(mocker):

    mock_make_request = mocker.patch("arctic3d.cli_localise.make_request")

    mock_make_request.return_value = {}

    result = call_uniprot("P12345")

    assert result == {}

    mock_make_request.side_effect = Exception("Mocked exception")

    result = call_uniprot("P12345")

    assert result is None
