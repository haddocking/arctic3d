from arctic3d.cli_localise import (
    main,
    get_quickgo_information,
    shorten_labels,
    get_uniprot_subcellular_location,
)

from pathlib import Path

import pytest
import os
import shutil

from . import golden_data


@pytest.fixture
def empty_cluster_filepath():
    """Empty cluster filepath."""
    return Path(golden_data, "clustered_interfaces_empty.out")


@pytest.fixture
def example_B_labels():
    """Example biological process labels."""
    return [
        "activation of cysteine-type endopeptidase activity involved in apoptotic process",  # noqa: E501
        "apoptotic signaling pathway",
        "positive regulation of transcription by RNA polymerase II",
    ]


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


def test_localise_cli_empty(empty_cluster_filepath):
    start_cwd = os.getcwd()
    run_dir = "arctic3d-localise"
    main(
        empty_cluster_filepath,
        "arctic3d-localise",
        None,
        None,
        None,
    )
    os.chdir(start_cwd)
    shutil.rmtree(run_dir)


def test_shorten_labels(example_B_labels):
    """Test shorten_labels."""
    obs_shortened_labels = shorten_labels(example_B_labels, 50)
    exp_shortened_labels = [
        "activation of cysteine-type endopeptidase activity...",
        "apoptotic signaling pathway",
        "positive regulation of transcription by RNA polymerase...",
    ]
    assert exp_shortened_labels == obs_shortened_labels


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
