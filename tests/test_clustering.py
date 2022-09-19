import os

import numpy as np
import pytest

from arctic3d.modules.clustering import (
    cluster_similarity_matrix,
    get_clustering_dict,
    get_residue_dict,
    write_clusters,
    write_residues,
)


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


def test_cluster_similarity_matrix():
    """Test correct clustering"""
    int_matrix = np.array([0.9972, 0.3742, 0.9736, 0.9996, 0.8841, 0.9991])
    ligands = ["int_1", "int_2", "int_3", "int_4"]
    clusters = cluster_similarity_matrix(int_matrix, ligands)
    expected_clusters = [1, 2, 1, 3]
    assert (clusters == expected_clusters).all()


def test_get_cl_dict():
    """Test correct retrieval of cl_dict."""
    clusters_list = [1, 1, 2, 3, 3, 4, 2]
    ligands_list = ["int1", "int2", "p53", "00", "int47", "antibody", "dimer"]
    expected_cl_dict = {
        1: ["int1", "int2"],
        2: ["p53", "dimer"],
        3: ["00", "int47"],
        4: ["antibody"],
    }
    observed_cl_dict = get_clustering_dict(clusters_list, ligands_list)
    assert expected_cl_dict, observed_cl_dict


def test_get_res_dict():
    """Test correct retrieval of res_dict."""
    interface_dict = {"int_1": [1, 2, 3], "int_2": [3, 4, 5], "int_3": [27, 28, 29]}
    cl_dict = {1: ["int_1", "int_2"], 2: ["int_3"]}
    expected_res_dict = {1: [1, 2, 3, 4, 5], 2: [27, 28, 29]}
    observed_res_dict = get_residue_dict(cl_dict, interface_dict)
    assert expected_res_dict == observed_res_dict


def test_write_clusters(reference_cl_dict):
    """Test write_clusters."""
    cl_filename = "clusters_test.out"
    write_clusters(reference_cl_dict, cl_filename)
    expected_content = (
        f"Cluster 1 -> int_1 int_2{os.linesep}Cluster 2 -> int_3{os.linesep}"
    )
    observed_content = open(cl_filename, "r").read()
    assert expected_content == observed_content
    os.unlink(cl_filename)


def test_write_residues(reference_res_dict):
    """Test write_residues."""
    res_filename = "residues_test.out"
    write_residues(reference_res_dict, res_filename)
    expected_content = (
        f"Cluster 1 -> 1 2 3 4 5{os.linesep}Cluster 2 -> 27 28 29{os.linesep}"
    )
    observed_content = open(res_filename, "r").read()
    assert expected_content == observed_content
    os.unlink(res_filename)
