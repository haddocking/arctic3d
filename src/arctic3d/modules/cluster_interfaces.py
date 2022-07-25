import logging

from arctic3d.modules.clustering import interface_clustering, write_residues
from arctic3d.modules.interface_matrix import interface_matrix

log = logging.getLogger("arctic3dlog")


def cluster_interfaces(interface_dict, pdb_path):
    """
    Wrapper to call interface_matrix and clustering

    Parameters
    ----------
    interface_dict : dict
        dictionary of all the interfaces (each one with its uniprot ID as key)
    pdb_path : str or Path
        pdb filename

    Returns
    -------
    clustered_residues : dict
        dictionary of the clustered interfaces
    """
    filtered_interfaces, matrix_path = interface_matrix(interface_dict, pdb_path)
    if len(filtered_interfaces) > 1:
        clustered_residues = interface_clustering(filtered_interfaces, matrix_path)
    elif len(filtered_interfaces) == 1:
        unique_interface = list(filtered_interfaces.keys())[0]
        cl_dict = {1: [unique_interface]}
        cl_filename = "clustered_residues.out"
        clustered_residues = write_residues(cl_dict, filtered_interfaces, cl_filename)
    else:
        clustered_residues = None
    return clustered_residues
