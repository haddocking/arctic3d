from arctic3d.modules.clustering import interface_clustering
from arctic3d.modules.interface_matrix import interface_matrix


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
    matrix_path = interface_matrix(interface_dict, pdb_path)
    clustered_residues = interface_clustering(matrix_path, interface_dict)
    return clustered_residues
