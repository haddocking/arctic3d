import MDAnalysis as mda
from scipy.spatial.distance import cdist
import numpy as np
import os

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
    print(f"computing scal_prod between {interface_one} and {interface_two}")
    len_one = len(interface_one)
    len_two = len(interface_two)
    #npairs = len_one * len_two
    #matrix_indices = np.zeros((npairs, 2 ), dtype=int)
    #matrix_idx = 0

    scalar_product = 0.0
    for res_one in range(len_one):
        for res_two in range(len_two):
            scalar_product += Jij_mat[interface_one[res_one], interface_two[res_two]]
            #idx = interface_one[res_one]*Jij_mat.shape[0] + interface_two[res_two]
            #print(f"idx {idx}")
            #matrix_indices[matrix_idx] = (res_one, res_two)
            #matrix_idx += 1
    print(f"scal prod {scalar_product}")
    #scalar_product = np.sum(Jij_mat[matrix_indices])
    return scalar_product

def get_coupling_matrix(mdu, int_resids):
    """
    Computes coupling matrix.

    Parameters
    ----------


    Returns
    -------
    """
    sel_residues = "name CA and resid " + " ".join([str(el) for el in int_resids])
    u = mdu.select_atoms(sel_residues)
    distmap = cdist(u.positions, u.positions)
    exp_factor = 4*SIGMA*SIGMA
    minus_dsq = - np.power(distmap,2)/(exp_factor)
    Jij_mat = np.exp(minus_dsq)
    #print(f"jij_mat {exp}")
    return Jij_mat

def output_interface_matrix(int_names, int_matrix, output_filename):
    """
    Writes the interface matrix
    Parameters
    ----------
    retained_interfaces : list
        list of the names of the interfaces
    
    int_matrix : 1D np.array
        interface distance matrix

    output_fl : str or Path

    Returns
    -------
    """
    matrix_idx = 0
    with open(output_filename, "w") as wmatrix:
        for int_one in range(len(int_names)):
            for int_two in range(int_one + 1, len(int_names)):
                string = f"{int_names[int_one]} {int_names[int_two]} {int_matrix[matrix_idx]:.4f}"
                string += os.linesep
                matrix_idx += 1
                wmatrix.write(string)

def get_unique_sorted_resids(interface_names, interface_dict):
    """

    """
    int_resids = []
    for key in interface_names:
        for resid in interface_dict[key]:
            if resid not in int_resids:
                int_resids.append(resid)
    int_resids.sort()
    return int_resids


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
    n_ret = len(retained_interfaces)
    print(f"{n_ret} retained_interfaces {retained_interfaces}")
    int_pairs = int(n_ret*(n_ret-1)/2)
    print(f"{int_pairs} pairs of interfaces")
    # getting all the residues
    int_resids = get_unique_sorted_resids(retained_interfaces, interface_dict)
    # using index to keep track of interfaces
    mapped_int_dict = {}
    for key in retained_interfaces:
        mapped_int_dict[key] = [int_resids.index(el) for el in interface_dict[key]]
    print(f"mapped_int_dict {mapped_int_dict}")
    # calculate coupling matrix
    Jij_mat = get_coupling_matrix(mdu, int_resids)
    # norms
    norms = np.zeros(n_ret)
    for key_idx in range(n_ret):
        key = retained_interfaces[key_idx]
        norms[key_idx] = compute_scalar_product(mapped_int_dict[key], mapped_int_dict[key], Jij_mat)
    #print(f"computed norms {norms}")
    # for each interface pair, calculate the scalar product
    scal_prods = np.zeros(int_pairs)
    prod_idx = 0
    for idx_one in range(n_ret):
        key_one = retained_interfaces[idx_one]
        for idx_two in range(idx_one + 1, n_ret):
            key_two = retained_interfaces[idx_two]
            scal_prods[prod_idx] = compute_scalar_product(mapped_int_dict[key_one],
                                                          mapped_int_dict[key_two],
                                                          Jij_mat
                                                         )
            prod_idx += 1
    #print(f"scalar products {scal_prods}")
    # calculate cosine and sine matrix
    cos_mat = np.zeros(int_pairs)
    mat_idx = 0
    for idx_one in range(n_ret):
        for idx_two in range(idx_one + 1, n_ret):
            cos_mat[mat_idx] = scal_prods[mat_idx]/np.sqrt(norms[idx_one]*norms[idx_two])
            mat_idx += 1
    sin_mat = np.ones(int_pairs) - np.power(cos_mat,2)
    out_fl = "interface.txt"
    output_interface_matrix(retained_interfaces, sin_mat, out_fl)

    
interface_dict = {#"int_1" : [1,2,5,6,9,10,11,14,15,16,18,19],
                  #"int_2" : [21,22,24,26,27,29,30,31,36,38,39,40],
                  "int_2" : [1,2],
                  "int_3" : [1,2,3,4],
                  "int_4" : [2,4]
                 }

interface_matrix(interface_dict, "/trinity/login/mgiulini/antibodies/4DN4-new/data/4DN4_l_u.pdb")