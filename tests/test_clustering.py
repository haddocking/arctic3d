import numpy as np

from arctic3d.modules.clustering import (  # write_clusters,; write_residues,
    cluster_distance_matrix,
)


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
