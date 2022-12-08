import logging

from arctic3d.modules.clustering import interface_clustering
from arctic3d.modules.interface_matrix import interface_matrix

log = logging.getLogger("arctic3d.log")


def cluster_interfaces(interface_dict, pdb_path, linkage_strategy, threshold):
    """
    Wrapper to call interface_matrix and clustering

    Parameters
    ----------
    interface_dict : dict
        dictionary of all the interfaces (each one with its uniprot ID as key)
    pdb_path : str or Path
        pdb filename
    linkage_strategy : str
        linkage strategy for clustering
    threshold : float
        threshold for clustering

    Returns
    -------
    clustered_residues : dict
        dictionary of the clustered interfaces
    """
    filtered_interfaces, matrix_path = interface_matrix(interface_dict, pdb_path)
    if len(filtered_interfaces) > 0:
        clustered_residues = interface_clustering(
            filtered_interfaces, matrix_path, linkage_strategy, threshold
        )
    else:
        clustered_residues = None
    return clustered_residues
