from pathlib import Path

import numpy as np
import pytest

from arctic3d.modules.clustering import (
    read_int_matrix,
    cluster_distance_matrix,
    write_clusters,
    write_residues,
)

from . import golden_data


def test_read_int_matrix_nonexisting():
    """Test error on non-existing path."""
    non_ex_path = "../dummy"
    with pytest.raises(Exception):
        read_int_matrix(non_ex_path)

def test_read_int_matrix():
    """Test correct reading of interface matrix."""
    matrix_path = Path(golden_data, "interface_matrix.txt")
    expected_int_matrix = np.array([0.9972, 0.3742, 0.9736, 0.9996, 0.8841, 0.9991])
    expected_ligands = ["int_1", "int_2", "int_3", "int_4"]
    observed_int_matrix, observed_ligands = read_int_matrix(matrix_path)
    assert expected_ligands == observed_ligands
    # now checking the matrix
    np.testing.assert_allclose(expected_int_matrix, observed_int_matrix, atol=0.0001)


def test_cluster_distance_matrix():
    """Test correct clustering"""
    int_matrix = np.array([0.9972, 0.3742, 0.9736, 0.9996, 0.8841, 0.9991])
    ligands = ["int_1", "int_2", "int_3", "int_4"]
    clusters = cluster_distance_matrix(int_matrix, ligands)
    expected_clusters = [1, 2, 1, 3]
    assert (clusters == expected_clusters).all()

def test_write_clusters():
    # TODO: implement this
    assert True
def test_write_residues():
    # TODO: implement this
    assert True