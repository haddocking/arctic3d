import logging
import os
import time

import MDAnalysis as mda
import numpy as np
from scipy.spatial.distance import cdist

SIGMA = 1.9
INTERFACE_COV_CUTOFF = 0.8

log = logging.getLogger("arctic3dlog")


def check_residues_coverage(interface, pdb_resids):
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
    filtered_interface = []
    for int_res in interface:
        if int_res in pdb_resids:
            coverage += 1.0
            filtered_interface.append(int_res)
    coverage = coverage / len(interface)
    return coverage, filtered_interface


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
    log.debug(f"computing scal_prod between {interface_one} and {interface_two}")
    len_one = len(interface_one)
    len_two = len(interface_two)
    scalar_product = 0.0
    for res_one in range(len_one):
        for res_two in range(len_two):
            scalar_product += Jij_mat[interface_one[res_one], interface_two[res_two]]
    return scalar_product


def get_coupling_matrix(mdu, int_resids):
    """
    Computes coupling matrix.

    Parameters
    ----------
    mdu : mdanalysis.universe
        universe of the reference pdb
    int_resids : list
        list of interacting residues

    Returns
    -------
    Jij_mat : np.array
        coupling matrix
    """
    sel_residues = "name CA and resid " + " ".join([str(el) for el in int_resids])
    u = mdu.select_atoms(sel_residues)
    if u.positions.shape[0] != len(int_resids):
        raise Exception(
            "shape mismatch: positions do not match input residues {int_resids}"
        )
    distmap = cdist(u.positions, u.positions)
    exp_factor = 4 * SIGMA * SIGMA
    minus_dsq = -np.power(distmap, 2) / (exp_factor)
    Jij_mat = np.exp(minus_dsq)
    return Jij_mat


def output_interface_matrix(int_names, int_matrix, output_filename):
    """
    Writes the interface matrix.

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
    log.info("Writing interface matrix on file {output_filename}")
    matrix_idx = 0
    with open(output_filename, "w") as wmatrix:
        for int_one in range(len(int_names)):
            for int_two in range(int_one + 1, len(int_names)):
                string = f"{int_names[int_one]} {int_names[int_two]} {int_matrix[matrix_idx]:.4f}"
                string += os.linesep
                matrix_idx += 1
                wmatrix.write(string)


def get_unique_sorted_resids(interface_dict):
    """
    Gets unique and sorted residues in an interface dictionary

    Parameters
    ----------
    interface_dict : dict
        dictionary of interfaces

    Returns
    -------
    int_resids : list
        list of sorted and unique residues
    """
    int_resids = []
    for key in interface_dict:
        for resid in interface_dict[key]:
            if resid not in int_resids:
                int_resids.append(resid)
    int_resids.sort()
    return int_resids


def filter_interfaces(interface_dict, pdb_resids):
    """
    Filters the interfaces accoriding to the residues present in the pdb

    Parameters
    ----------
    interface_dict : dict
        dictionary of interfaces

    pdb_resids : np.array
        residues present in the pdb

    Returns
    -------
    retained_interfaces : dict
        dictionary of the retained and filtered interfaces
        example : interface_dict = {"a" : [1,2], "b" : [2,3,4], "c": [5,6,7]}
                  pdb_resids = np.array([3,4,5,6,7])
        then, if INTERFACE_COV_CUTOFF < 0.66:
            retained_interfaces = {"b": [3,4], "c" : [5,6,7]}
        else:
            retained_interfaces = {"c" : [5,6,7]}
    """
    log.debug("filtering interface dictionary")
    retained_interfaces = {}
    for key in interface_dict.keys():
        coverage, filtered_interface = check_residues_coverage(
            interface_dict[key], pdb_resids
        )
        if coverage > INTERFACE_COV_CUTOFF:
            retained_interfaces[key] = filtered_interface
    log.info(f"{len(retained_interfaces.keys())} retained_interfaces")
    return retained_interfaces


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
    start_time = time.time()
    log.info("computing interface matrix")
    log.info(f"{os.listdir()}")
    if not os.path.exists(pdb_path):
        raise Exception(f"pdb_path {pdb_path} does not exist")
    mdu = mda.Universe(pdb_path)
    pdb_resids = mdu.select_atoms("name CA").resids
    retained_interfaces = filter_interfaces(interface_dict, pdb_resids)
    ret_keys = list(retained_interfaces.keys())
    n_ret = len(ret_keys)
    int_pairs = int(n_ret * (n_ret - 1) / 2)
    log.info(f"{int_pairs} pairs of interfaces")
    # getting all the residues
    int_resids = get_unique_sorted_resids(retained_interfaces)
    log.info(f"interacting residues {int_resids}")
    # using index to keep track of interfaces
    mapped_int_dict = {}
    for key in ret_keys:
        mapped_int_dict[key] = [int_resids.index(el) for el in retained_interfaces[key]]
    # calculate coupling matrix
    Jij_mat = get_coupling_matrix(mdu, int_resids)
    # norms
    norms = np.zeros(n_ret)
    for key_idx in range(n_ret):
        norms[key_idx] = compute_scalar_product(
            mapped_int_dict[ret_keys[key_idx]],
            mapped_int_dict[ret_keys[key_idx]],
            Jij_mat,
        )
    # for each interface pair, calculate the scalar product
    scal_prods = np.zeros(int_pairs)
    prod_idx = 0
    for idx_one in range(n_ret):
        for idx_two in range(idx_one + 1, n_ret):
            scal_prods[prod_idx] = compute_scalar_product(
                mapped_int_dict[ret_keys[idx_one]],
                mapped_int_dict[ret_keys[idx_two]],
                Jij_mat,
            )
            prod_idx += 1
    # calculate cosine and sine matrix
    cos_mat = np.zeros(int_pairs)
    mat_idx = 0
    for idx_one in range(n_ret):
        for idx_two in range(idx_one + 1, n_ret):
            cos_mat[mat_idx] = scal_prods[mat_idx] / np.sqrt(
                norms[idx_one] * norms[idx_two]
            )
            mat_idx += 1
    sin_mat = np.ones(int_pairs) - np.power(cos_mat, 2)
    out_fl = "interface_matrix.txt"
    output_interface_matrix(ret_keys, sin_mat, out_fl)
    elap_time = round((time.time() - start_time), 2)
    log.info(f"interface matrix calculated in {elap_time} seconds")
    return out_fl
