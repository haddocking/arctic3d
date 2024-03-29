import os
from pathlib import Path

import pytest

from arctic3d.modules.output import (
    create_output_folder,
    output_pdb,
    read_residues_probs,
    remove_duplicate_labels,
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
def clust_probs_example():
    """Example file with residue probabilities."""
    return Path(golden_data, "clustered_residues_probs.out")


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
    res_probs = {1: {3: 0.2, 4: 0.75}, 2: {27: 1.0, 28: 1.0}}
    full_resnames_dict = {3: "ALA", 4: "GLY", 27: "ARG", 28: "LYS"}
    expected_content = (
        f"Cluster 1 : 2 residues{os.linesep}"
        f"rank\tresid\tresname\tprobability{os.linesep}"
        f"1\t4\tGLY\t0.750{os.linesep}"
        f"2\t3\tALA\t0.200{os.linesep}"
        f"{os.linesep}"
        f"Cluster 2 : 2 residues{os.linesep}"
        f"rank\tresid\tresname\tprobability{os.linesep}"
        f"1\t27\tARG\t1.000{os.linesep}"
        f"2\t28\tLYS\t1.000{os.linesep}"
        f"{os.linesep}"
    )
    res_probs_filename = "residues_probs_test.out"
    write_residues_probs(res_probs, res_probs_filename, full_resnames_dict)
    observed_content = open(res_probs_filename, "r").read()
    assert expected_content == observed_content
    os.unlink(res_probs_filename)


def test_read_res_probs(clust_probs_example):
    """Test read_residues_probs."""
    example_res_probs = {
        2: {
            36: 1.000,
            37: 1.000,
            95: 1.000,
            99: 1.000,
            211: 1.000,
            233: 1.000,
            238: 1.000,
            237: 0.857,
            464: 0.714,
            241: 0.571,
            365: 0.571,
            375: 0.571,
            32: 0.429,
            33: 0.429,
            96: 0.286,
            363: 0.286,
        },
        1: {47: 1.0, 48: 1.0, 51: 1.00},
    }
    observed_res_probs = read_residues_probs(clust_probs_example)
    assert example_res_probs == observed_res_probs


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


def test_create_output_folder():
    """Test if the expected run_dir is effectively created."""
    uniprot_id = "fake_uniprot"
    create_output_folder(output_dir=None, uniprot_id=uniprot_id)
    exp_run_dir = Path(f"arctic3d-{uniprot_id}")
    assert Path.exists(exp_run_dir)
    os.rmdir(exp_run_dir)


def test_setup_output_folder(inp_pdb):
    """Test the correct setup of the output folder."""
    run_dir = "dummy_output"
    start_cwd = os.getcwd()
    create_output_folder(run_dir)
    input_files = {"pdb": inp_pdb}
    setup_output_folder(run_dir, input_files)
    obs_cwd = Path(os.getcwd())
    exp_cwd = Path(start_cwd, run_dir)
    assert exp_cwd == obs_cwd
    os.chdir(start_cwd)
    assert Path.exists(Path(run_dir, "input_data"))
    assert Path.exists(Path(run_dir, "input_data", inp_pdb.name))
    os.unlink(Path(run_dir, "input_data", inp_pdb.name))
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


def test_remove_duplicate_labels():
    """Test remove_duplicate_labels."""
    tmp_labels = ["Polymerase...", "Polymerase...", "Polymerase..."]
    tmp_values = [2, 3, 1]
    exp_labels = ["Polymerase..."]
    exp_values = [2]
    obs_labels, obs_values = remove_duplicate_labels(tmp_labels, tmp_values)
    assert exp_labels == obs_labels
    assert exp_values == obs_values
