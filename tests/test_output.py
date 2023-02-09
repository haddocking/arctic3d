import os
from pathlib import Path

import pytest

from arctic3d.modules.output import (
    output_pdb,
    setup_output_folder,
    shorten_labels,
    write_dict,
    write_residues_probs,
)

from . import golden_data


@pytest.fixture
def inp_pdb():
    return Path(golden_data, "1rypB_r_b.pdb")


@pytest.fixture
def reference_cl_dict():
    """Reference dictionary of clustered interfaces."""
    cl_dict = {1: ["int_1", "int_2"], 2: ["int_3"]}
    return cl_dict


@pytest.fixture
def reference_res_dict():
    """Reference dictionary of clustered interfaces."""
    res_dict = {1: [1, 2, 3, 4, 5], 2: [27, 28, 29]}
    return res_dict


@pytest.fixture
def example_B_labels():
    """Example biological process labels."""
    return [
        "activation of cysteine-type endopeptidase activity involved in apoptotic process",  # noqa: E501
        "apoptotic signaling pathway",
        "positive regulation of transcription by RNA polymerase II",
    ]


def test_write_clusters(reference_cl_dict):
    """Test write_dict for clusters."""
    cl_filename = "clusters_test.out"
    write_dict(reference_cl_dict, cl_filename, keyword="Cluster")
    expected_content = (
        f"Cluster 1 -> int_1 int_2{os.linesep}Cluster 2 -> int_3{os.linesep}"
    )
    observed_content = open(cl_filename, "r").read()
    assert expected_content == observed_content
    os.unlink(cl_filename)


def test_write_residues(reference_res_dict):
    """Test write_dict for residues."""
    res_filename = "residues_test.out"
    write_dict(reference_res_dict, res_filename, keyword="Cluster")
    expected_content = (
        f"Cluster 1 -> 1 2 3 4 5{os.linesep}Cluster 2 -> 27 28 29{os.linesep}"
    )
    observed_content = open(res_filename, "r").read()
    assert expected_content == observed_content
    os.unlink(res_filename)


def test_write_interfaes(reference_res_dict):
    """Test write_dict for interfaces."""
    res_filename = "residues_test.out"
    write_dict(reference_res_dict, res_filename, keyword="Interface")
    expected_content = (
        f"Interface 1 -> 1 2 3 4 5{os.linesep}Interface 2 -> 27 28"
        f" 29{os.linesep}"
    )
    observed_content = open(res_filename, "r").read()
    assert expected_content == observed_content
    os.unlink(res_filename)


def test_write_res_probs():
    """Test write_residue_probs."""
    example_res_probs = {1: {3: 0.2, 4: 0.75}, 2: {27: 1.0, 28: 1.0}}
    expected_content = (
        f"Cluster 1 : 2 residues{os.linesep}"
        f"rank\tresid\tprobability{os.linesep}"
        f"1\t4\t0.750{os.linesep}"
        f"2\t3\t0.200{os.linesep}"
        f"{os.linesep}"
        f"Cluster 2 : 2 residues{os.linesep}"
        f"rank\tresid\tprobability{os.linesep}"
        f"1\t27\t1.000{os.linesep}"
        f"2\t28\t1.000{os.linesep}"
        f"{os.linesep}"
    )
    res_probs_filename = "residues_probs_test.out"
    write_residues_probs(example_res_probs, res_probs_filename)
    observed_content = open(res_probs_filename, "r").read()
    assert expected_content == observed_content
    os.unlink(res_probs_filename)


def test_output_pdb(inp_pdb):
    """Test output_pdb."""
    example_res_probs = {1: {3: 0.2, 4: 0.75}}
    output_files = output_pdb(inp_pdb, example_res_probs)
    original_content = open(inp_pdb, "r").read().split(os.linesep)
    original_content = [el for el in original_content if el.startswith("ATOM")]
    # check file existence
    assert Path.exists(output_files[0])
    observed_content = open(output_files[0], "r").read().split(os.linesep)
    observed_content = [el for el in observed_content if el.startswith("ATOM")]
    # assert equal length
    assert len(original_content) == len(observed_content)
    # assert equal content (except for b factors)
    for ln_id in range(len(observed_content)):
        assert original_content[ln_id][:60] == observed_content[ln_id][:60]
    os.unlink(output_files[0])


def test_run_dir():
    """Test if the expected run_dir is effectively created."""
    run_dir = "run_dir"
    uniprot_id = "fake_uniprot"
    start_cwd = os.getcwd()
    setup_output_folder(uniprot_id, [], run_dir)
    obs_cwd = Path(os.getcwd())
    exp_cwd = Path(start_cwd, run_dir)
    assert exp_cwd == obs_cwd
    os.chdir(start_cwd)
    os.rmdir(Path(run_dir, "input_data"))
    os.rmdir(run_dir)


def test_shorten_labels(example_B_labels):
    """Test shorten_labels."""
    obs_shortened_labels = shorten_labels(example_B_labels, 50)
    exp_shortened_labels = [
        "activation of cysteine-type endopeptidase activity...",
        "apoptotic signaling pathway",
        "positive regulation of transcription by RNA polymerase...",
    ]
    assert exp_shortened_labels == obs_shortened_labels
