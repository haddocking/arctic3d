import os
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pytest

from arctic3d.modules.interface_matrix import (
    check_residues_coverage,
    compute_scalar_product,
    filter_interfaces,
    get_coupling_matrix,
    interface_matrix,
    read_int_matrix,
)

from . import golden_data


@pytest.fixture
def example_pdbpath():
    """Example pdb path."""
    return Path(golden_data, "1rypB_r_b.pdb")


@pytest.fixture
def example_mdu(example_pdbpath):
    """Example mdanalysis universe."""
    return mda.Universe(example_pdbpath)


@pytest.fixture
def example_interface_dict():
    """Example interface dictionary."""
    interface_dict = {"int_1": [1, 2], "int_2": [1, 2, 4], "int_3": [250, 251]}
    return interface_dict


@pytest.fixture
def reference_jij():
    jij = np.array(
        [
            [1.0, 0.358133, 0.031553],
            [0.358133, 1.0, 0.366509],
            [0.031553, 0.366509, 1.0],
        ]
    )
    return jij


def test_check_residues_coverage():
    """Test check_residues_coverage."""
    interface_one = [1, 2, 3]
    interface_two = [2, 3, 4, 5]
    pdb_resids = [1, 2, 3, 4]
    cov_one, filtered_int_one = check_residues_coverage(
        interface_one, pdb_resids
    )
    assert cov_one == 1.0
    assert filtered_int_one == interface_one

    cov_two, filtered_int_two = check_residues_coverage(
        interface_two, pdb_resids
    )
    assert cov_two == 0.75
    expected_filtered_int_two = [2, 3, 4]
    assert expected_filtered_int_two == filtered_int_two


def test_get_coupling_matrix_empty(example_mdu):
    """Test get_coupling_matrix with a single residue"""
    int_resids = [1]
    observed_jij = get_coupling_matrix(example_mdu, int_resids)
    expected_jij = np.array([1.0])
    assert observed_jij == expected_jij


def test_get_coupling_matrix(example_mdu, reference_jij):
    """Test get_coupling_matrix with a set of residues"""
    int_resids = [1, 2, 3]
    observed_jij = get_coupling_matrix(example_mdu, int_resids)
    np.testing.assert_allclose(reference_jij, observed_jij, atol=0.00001)


def test_compute_scalar_product(reference_jij):
    """Test compute_scalar_product."""
    interface_one = [0, 1, 2]
    observed_norm = compute_scalar_product(
        interface_one, interface_one, reference_jij
    )
    expected_norm = 4.51239
    np.testing.assert_allclose(expected_norm, observed_norm, atol=0.00001)
    interface_two = [0, 1]
    observed_scal_prod = compute_scalar_product(
        interface_one, interface_two, reference_jij
    )
    expected_scal_prod = 3.11433
    np.testing.assert_allclose(
        expected_scal_prod, observed_scal_prod, atol=0.00001
    )


def test_filter_interfaces(example_mdu, example_interface_dict):
    """Test filter_interfaces."""
    expected_filter_dict = {"int_1": [1, 2], "int_2": [1, 2, 4]}
    pdb_resids = example_mdu.select_atoms("name CA").resids
    observed_filter_dict = filter_interfaces(
        example_interface_dict, pdb_resids
    )
    assert expected_filter_dict == observed_filter_dict


def test_interface_matrix(example_interface_dict, example_pdbpath):
    """Test interface_matrix"""
    # defining expected quantities
    expected_interface_matrix = np.array([0.2515])
    expected_int_filename = "interface_matrix.txt"
    expected_filter_dict = {"int_1": [1, 2], "int_2": [1, 2, 4]}
    expected_interface_names = ["int_1", "int_2"]
    # calculate interface matrix
    observed_filter_ints, observed_int_filename = interface_matrix(
        example_interface_dict, example_pdbpath
    )
    assert expected_filter_dict == observed_filter_ints
    assert expected_int_filename == observed_int_filename
    # read interface matrix
    observed_int_matrix, obs_int_names = read_int_matrix(observed_int_filename)
    np.testing.assert_allclose(
        expected_interface_matrix, observed_int_matrix, atol=0.00001
    )
    assert expected_interface_names == obs_int_names
    os.unlink(observed_int_filename)


def test_read_int_matrix_nonexisting():
    """Test error on non-existing path."""
    non_ex_path = "../dummy"
    with pytest.raises(Exception):
        read_int_matrix(non_ex_path)


def test_read_int_matrix():
    """Test correct reading of interface matrix."""
    matrix_path = Path(golden_data, "interface_matrix.txt")
    expected_int_matrix = np.array(
        [0.9972, 0.3742, 0.9736, 0.9996, 0.8841, 0.9991]
    )
    expected_ligands = ["int_1", "int_2", "int_3", "int_4"]
    observed_int_matrix, observed_ligands = read_int_matrix(matrix_path)
    assert expected_ligands == observed_ligands
    # now checking the matrix
    np.testing.assert_allclose(
        expected_int_matrix, observed_int_matrix, atol=0.0001
    )
