import MDAnalysis as mda
from scipy.spatial.distance import pdist
import numpy as np

SIGMA = 1.9
INTERFACE_COV_CUTOFF = 0.8

def check_residues(interface, pdb_resids):
    """
    Checks if the residues of the interface are included in the pdb_resids list

    It does not do any numbering check!

    Parameters
    ----------
    interface : list
        list of interface residues
    pdb_residues : list or np.array
        list of all residues contained in a pdb

    Returns
    -------
    coverage : float
        fraction of resids of the interface found in the pdb
    """
    coverage = 0.0
    for int_res in interface:
        if int_res in pdb_resids:
            coverage += 1.0
    coverage = coverage/len(interface)
    return coverage

def compute_norm(interface, Jij_mat):
    """
    Computes the norm of an interface.

    Parameters
    ----------
    interface : list
        list of interface residues
    Jij_mat : np.array
        coupling matrix

    Returns
    -------
    norm : float
        interface norm
    """
    len_int = len(interface)
    norm = float(len_int) # self-interaction term
    knt = 0
    npairs = len_int * (len_int - 1)/2
    for i in range(npairs):
        norm += Jij_mat[i]
    return norm

def compute_scalar_product(interface_one, interface_two, Jij_mat):
    """
    Computes the scalar product between two interfaces.

    Parameters
    ----------
    interface_one : list
        list of interface residues
    interface_two : list
        list of interface residues
    Jij_mat : np.array
        coupling matrix

    Returns
    -------
    scalar_product : float
        scalar product between the two interfaces
    """
    len_one = len(interface_one)
    len_two = len(interface_two)
    npairs = len_one * len_two
    matrix_indices = np.array(npairs)
    for res_one in range(len_one):
        for res_two in range(len_two):
            idx = interface_one[res_one]*Jij_mat.shape[0] + interface_one[res_two]
            matrix_indices[matrix_idx] = idx
            matrix_idx += 1
    print(f"matrix_indices {matrix_indices}")
    scalar_product = np.sum(Jij_mat[matrix_indices])
    return scalar_product


def interface_matrix(interface_dict, pdb_path):
    """
    Computes the interface matrix.

    Parameters
    ----------
    interface_dict : dict
        dictionary of all the interfaces (each one with its uniprot ID as key)
    pdb_path : str or Path
        path to the pdb of interest

    Returns
    -------
    interface_matrix : np.array
        interface matrix
    """
    mdu = mda.Universe(pdb_path)
    pdb_resids = mdu.select_atoms('name CA').resids
    retained_interfaces = [] # list of uniprot IDs whose interface is part of the ref pdb
    for key in interface_dict.keys():
        chk = check_residues(interface_dict[key], pdb_resids)
        if chk:
            retained_interfaces.append(key)
    
    int_resids = []
    for key in retained_interfaces:
        for resid in interface_dict[key]:
            if resid not in int_resids:
                int_resids.append(resid)
    int_resids.sort()
    
    sel_residues = "name CA and resid " + " ".join([str(el) for el in int_resids])
    print(sel_residues)
    u = mdu.select_atoms(sel_residues)
    distmap = pdist(u.positions)

    exp_factor = 4*1.9*1.9
    minus_dsq = - np.power(distmap,2)/(exp_factor)
    exp = np.exp(minus_dsq)

    mapped_int_dict = {}
    for key in retained_interfaces:
        mapped_int_dict[key] = [interface_dict[key].index(el) for el in interface_dict[key]]
    
    norms = np.zeros(len(retained_interfaces))
    for key_idx in range(len(retained_interfaces)):
        norms[key_idx] = compute_norm(mapped_int_dict[key_idx], exp)
    print(f"computed norms {norms}")
    
    # for each pair, calculate the scalar product
    scal_prods = np.zeros(len(retained_interfaces))
    
