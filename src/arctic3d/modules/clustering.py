"""Clustering module."""

import logging
import os
import time

import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage

from arctic3d.modules.interface_matrix import read_int_matrix

LINKAGE = "average"
# THRESHOLD = 0.7071  # np.sqrt(2)/2
THRESHOLD = 0.8660  # np.sqrt(3)/2

log = logging.getLogger("arctic3dlog")


def cluster_similarity_matrix(int_matrix, entries, plot=False):
    """
    Does the clustering.

    Parameters
    ----------
    int_matrix : np.array
        1D condensed interface similarity matrix
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
        dendrogram(Z, color_threshold=THRESHOLD, labels=entries)
        plt.xlabel("Interface Names")
        plt.ylabel("Similarity")
        plt.savefig(dendrogram_figure_filename)
        plt.close()
    # clustering
    clusters = fcluster(Z, t=THRESHOLD, criterion="distance")
    log.info("Dendrogram created and clustered.")
    log.debug(f"Clusters = {clusters}")
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
    log.info(f"Cluster dictionary {cl_dict}")
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


def interface_clustering(interface_dict, matrix_filename):
    """
    Clusters the interface matrix.

    Parameters
    ----------
    interface_dict : dict
        dictionary of all the interfaces (each one with its uniprot ID as key)

    matrix_filename : str or Path
        input interface matrix

    Returns
    -------
    clustered_residues : dict
    """
    start_time = time.time()
    log.info("Clustering interface matrix")
    # check if there's only a single interface
    if len(interface_dict) == 1:
        clusters = [1]
        entries = list(interface_dict.keys())  # the only entry
    else:
        int_matrix, entries = read_int_matrix(matrix_filename)  # read matrix
        # cluster matrix. TODO: make plot an external parameter
        clusters = cluster_similarity_matrix(int_matrix, entries, plot=True)
    # write clusters
    cl_filename = "clustered_interfaces.out"
    cl_dict = write_clusters(clusters, entries, cl_filename)
    # write clustered residues
    res_filename = "clustered_residues.out"
    clustered_residues = write_residues(cl_dict, interface_dict, res_filename)
    # write time
    elap_time = round((time.time() - start_time), 3)
    log.info(f"Clustering performed in {elap_time} seconds")
    return clustered_residues
