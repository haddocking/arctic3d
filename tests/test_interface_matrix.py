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
)

from . import golden_data


@pytest.fixture
def example_mdu():
    """Example mdanalysis universe."""
    pdb_path = Path(golden_data, "1rypB_r_b.pdb")
    return mda.Universe(pdb_path)


@pytest.fixture
def example_interface_dict():
    """Example interface dictionary."""
    interface_dict = {"int_1": [1, 2], "int_2": [1, 2, 4], "int_3": [250, 251]}
    return interface_dict


def test_check_residues_coverage():
    """Test check_residues_coverage."""
    interface_one = [1, 2, 3]
    interface_two = [2, 3, 4, 5]
    pdb_resids = [1, 2, 3, 4]
    cov_one, filtered_int_one = check_residues_coverage(interface_one, pdb_resids)
    assert cov_one == 1.0
    assert filtered_int_one == interface_one

    cov_two, filtered_int_two = check_residues_coverage(interface_two, pdb_resids)

    assert cov_two == 0.75
    expected_filtered_int_two = [2, 3, 4]
    assert expected_filtered_int_two == filtered_int_two


def test_get_coupling_matrix_empty(example_mdu):
    """Test get_coupling_matrix with a single residue"""
    int_resids = [1]
    observed_jij = get_coupling_matrix(example_mdu, int_resids)
    expected_jij = np.array([1.0])
    assert observed_jij == expected_jij


def test_get_coupling_matrix(example_mdu):
    """Test get_coupling_matrix with a set of residues"""
    int_resids = [1, 2, 3]
    expected_jij = np.array(
        [
            [1.0, 0.358133, 0.031553],
            [0.358133, 1.0, 0.366509],
            [0.031553, 0.366509, 1.0],
        ]
    )
    observed_jij = get_coupling_matrix(example_mdu, int_resids)
    np.testing.assert_allclose(expected_jij, observed_jij, atol=0.00001)


def test_compute_scalar_product():
    """Test compute_scalar_product."""
    jij = np.array(
        [
            [1.0, 0.358133, 0.031553],
            [0.358133, 1.0, 0.366509],
            [0.031553, 0.366509, 1.0],
        ]
    )
    interface_one = [0, 1, 2]
    observed_norm = compute_scalar_product(interface_one, interface_one, jij)
    expected_norm = 4.51239
    np.testing.assert_allclose(expected_norm, observed_norm, atol=0.00001)
    interface_two = [0, 1]
    observed_scal_prod = compute_scalar_product(interface_one, interface_two, jij)
    expected_scal_prod = 3.11433
    np.testing.assert_allclose(expected_scal_prod, observed_scal_prod, atol=0.00001)


def test_filter_interfaces(example_mdu, example_interface_dict):
    """Test filter_interfaces."""
    expected_filter_dict = {"int_1": [1, 2], "int_2": [1, 2, 4]}
    pdb_resids = example_mdu.select_atoms("name CA").resids
    observed_filter_dict = filter_interfaces(example_interface_dict, pdb_resids)
    assert expected_filter_dict == observed_filter_dict
