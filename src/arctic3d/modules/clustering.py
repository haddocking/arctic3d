"""Clustering module."""

import logging
import numpy as np
import pandas as pd
import os
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram
import matplotlib.pyplot as plt
import time


LINKAGE = "single"
THRESHOLD = 0.7071 # np.sqrt(2)/2

log = logging.getLogger("arctic3dlog")

def read_int_matrix(filename):
    """
    Read the interface matrix.
    
    Parameters
    ----------
    filename : str or Path
        interface matrix filename
    Returns
    -------
    int_matrix : np.array
        interface matrix
    """
    int_matrix = pd.read_csv(filename, header=None, sep=" ")
    int_matrix.columns = ["lig1", "lig2", "D"]
    # first check: it must be a 1D condensed distance matrix
    nligands = 0.5 + np.sqrt(0.25 + 2*int_matrix.shape[0])
    int_nligands = int(nligands)
    if abs(nligands - int_nligands) > 0.00001:
        raise Exception(f"npairs {int_matrix.shape[0]}: interface matrix should be a 1D condensed distance matrix")
    ligand_names = [int_matrix.iloc[0,0]]
    for lig in np.unique(int_matrix.iloc[:,1]):
        ligand_names.append(lig)
    log.info(f"ordered ligand names {ligand_names}")
    return int_matrix.iloc[:,2], ligand_names

def cluster_distance_matrix(int_matrix, entries, plot=False):
    """
    Does the clustering.
    
    Parameters
    ----------
    int_matrix : np.array
        1D condensed interface matrix
    entries : list
        names of the ligands
    plot : bool
        if True, plot the dendrogram
    Returns
    -------
    clusters : list
        list of clusters ID, each one associated to an entry
    """
    Z = linkage(int_matrix, LINKAGE)
    if plot:
        dendrogram_figure_filename = "dendrogram_" + LINKAGE + ".png"
        plt.figure()
        dn = dendrogram(
            Z,
            color_threshold=THRESHOLD,
            labels=entries
        )
        plt.savefig(dendrogram_figure_filename)
        plt.close()
    # clustering
    clusters = fcluster(Z, t = THRESHOLD, criterion="distance")
    log.info(f"dendrogram created and clustered. Clusters = {clusters}")
    return clusters

def write_clusters(clusters, ligands, cl_filename):
    """
    Writes clusters to file.
    
    Parameters
    ----------
    clusters : list
        list of cluster IDs
    ligands : list
        names of the ligands
    cl_filename : str or Path
        name of the output filename

    Returns
    -------
    cl_dict : dict
        dictionary of clustered interfaces
        example { 1 : ['interface_1', 'interface_3'] ,
                  2 : ['interface_2'],
                  ...
                }
    """
    log.info(f"Writing clusters to file {cl_filename}")
    cl_dict = {}
    #
    for cl in range(len(clusters)):
        if clusters[cl] not in cl_dict.keys():
            cl_dict[clusters[cl]] = [ligands[cl]]
        else:
            cl_dict[clusters[cl]].append(ligands[cl])
    with open(cl_filename, "w") as wfile:
        for key in cl_dict.keys():
            cl_string = " ".join(cl_dict[key])
            wfile.write(f"Cluster {key} -> " + cl_string + os.linesep)
    return cl_dict


def write_residues(cl_dict, interface_dict, res_filename):
    """
    Writes the clustered residues to file

    cl_dict : dict
        dictionary of the clustered interfaces
    interface_dict : dict
        dictionary of all the interfaces (each one with its uniprot ID as key)
    res_filename : str or Path
        output filename 
    Returns
    -------
    cl_residues : dict
        dictionary of clustered residues
    """
    clustered_residues = {}
    for key in cl_dict.keys():
        residues = []
        for int_id in cl_dict[key]:
            residues.extend(interface_dict[int_id])
        unique_cl_residues = list(set(residues))
        unique_cl_residues.sort()
        clustered_residues[key] = unique_cl_residues
    # write to file
    with open(res_filename, "w") as wfile:
        for key in clustered_residues.keys():
            cl_string = " ".join([str(el) for el in clustered_residues[key]])
            wfile.write(f"Cluster {key} -> " + cl_string + os.linesep)
    return clustered_residues

def interface_clustering(matrix_filename, interface_dict):
    """
    Clusters the interface matrix.
    
    Parameters
    ----------
    matrix_filename : str or Path
        input interface matrix
    interface_dict : dict
        dictionary of all the interfaces (each one with its uniprot ID as key)
    """
    start_time = time.time()
    log.info("clustering interface matrix")
    # read matrix
    if os.path.exists(matrix_filename):
        int_matrix, entries = read_int_matrix(matrix_filename)
    else:
        raise Exception(f"input path {matrix_filename} does not exist!")
    # cluster matrix
    # TODO: make plot an external parameter
    clusters = cluster_distance_matrix(int_matrix, entries, plot=True)
    # write clusters
    cl_filename = "clustered_interfaces.out"
    cl_dict = write_clusters(clusters, entries, cl_filename)
    # write clustered residues
    res_filename = "clustered_residues.out"
    write_residues(cl_dict, interface_dict, res_filename)
    # write time
    elap_time = round((time.time() - start_time), 3)
    log.info(f"clustering performed in {elap_time} seconds")
