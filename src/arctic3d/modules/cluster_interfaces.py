import logging

from arctic3d.modules.clustering import interface_clustering
from arctic3d.modules.interface_matrix import interface_matrix

log = logging.getLogger("arctic3dlog")


def cluster_interfaces(interface_dict, pdb_path, filter=True):
    """
    Wrapper to call interface_matrix and clustering

    Parameters
    ----------
    interface_dict : dict
        dictionary of all the interfaces (each one with its uniprot ID as key)
    pdb_path : str or Path
        pdb filename
    filter : bool
        filter the interfaces. Could be already filtered.

    Returns
    -------
    clustered_residues : dict
        dictionary of the clustered interfaces
    """
    filtered_interfaces, matrix_path = interface_matrix(
        interface_dict, pdb_path, filter
    )
    if len(filtered_interfaces) > 0:
        clustered_residues = interface_clustering(filtered_interfaces, matrix_path)
    else:
        clustered_residues = None
    return clustered_residues
