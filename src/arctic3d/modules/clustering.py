"""Clustering module."""

import logging
import os
import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage

from arctic3d.modules.interface_matrix import read_int_matrix

LINKAGE = "average"
# THRESHOLD = 0.7071  # np.sqrt(2)/2
THRESHOLD = 0.8660  # np.sqrt(3)/2

log = logging.getLogger("arctic3dlog")


def cluster_similarity_matrix(int_matrix, entries, threshold=THRESHOLD, plot=False):
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
    log.info(f"Clustering with threshold {threshold}")
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
    clusters = fcluster(Z, t=threshold, criterion="distance")
    log.info("Dendrogram created and clustered.")
    log.debug(f"Clusters = {clusters}")
    return clusters


def get_clustering_dict(clusters, ligands):
    """
    Gets dictionary of clusters.

    Parameters
    ----------
    clusters : list
        list of cluster IDs
    ligands : list
        names of the ligands

    Returns
    -------
    cl_dict : dict
        dictionary of clustered interfaces
        example { 1 : ['interface_1', 'interface_3'] ,
                  2 : ['interface_2'],
                  ...
                }
    """
    cl_dict = {}
    # loop over clusters
    for cl in range(len(clusters)):
        if clusters[cl] not in cl_dict.keys():
            cl_dict[clusters[cl]] = [ligands[cl]]
        else:
            cl_dict[clusters[cl]].append(ligands[cl])
    log.info(f"Cluster dictionary {cl_dict}")
    return cl_dict


def write_clusters(cl_dict, cl_filename):
    """
    Writes clusters to file.

    Parameters
    ----------
    cl_dict : dict
        dictionary of clustered interfaces
    cl_filename : str or Path
        name of the output filename
    """
    log.info(f"Writing clusters to file {cl_filename}")
    with open(cl_filename, "w") as wfile:
        for key in cl_dict.keys():
            cl_string = " ".join(cl_dict[key])
            wfile.write(f"Cluster {key} -> " + cl_string + os.linesep)


def get_residue_dict(cl_dict, interface_dict):
    """
    Gets dictionary of clustered residues.

    Parameters
    ----------
    cl_dict : dict
        dictionary of the clustered interfaces
    interface_dict : dict
        dictionary of all the interfaces (each one with its uniprot ID as key)

    Returns
    -------
    clustered_residues : dict
        dictionary of clustered residues
        example { 1 : [1,2,3,5,6,8] ,
                  2 : [29,30,31],
                  ...
                }
    cl_residues_probs : dict of dicts
        dictionary of probabilities for clustered residues
        example { 1 : {1:0.7, 2:0.2, 3:0.4 ...}
                  ...
                }
    """
    clustered_residues = {}
    cl_residues_probs = {}
    for key in cl_dict.keys():
        denom = len(cl_dict[key])
        residues = []
        for int_id in cl_dict[key]:
            residues.extend(interface_dict[int_id])
        unique_res = np.unique(residues, return_counts=True)
        cl_residues_probs[key] = {}
        # assign probabilities
        for res_idx, res in enumerate(unique_res[0]):
            res_prob = unique_res[1][res_idx] / denom
            cl_residues_probs[key][res] = res_prob
        clustered_residues[key] = list(unique_res[0])
    return clustered_residues, cl_residues_probs


def write_residues(res_dict, res_filename):
    """
    Writes clustered residues to file.

    Parameters
    ----------
    res_dict : dict
        dictionary of clustered residues
    res_filename : str or Path
        output filename
    """
    # write to file
    with open(res_filename, "w") as wfile:
        for key in res_dict.keys():
            cl_string = " ".join([str(el) for el in res_dict[key]])
            wfile.write(f"Cluster {key} -> " + cl_string + os.linesep)


def write_residues_probs(cl_residues_probs, res_probs_filename):
    """
    Writes clustered residues to file with their probability.

    Parameters
    ----------
    res_dict : dict
        dictionary of clustered residues
    res_probs_filename : str or Path
        output filename
    """
    # write to file
    with open(res_probs_filename, "w") as wfile:
        for key in cl_residues_probs.keys():
            cl_string = f"Cluster {key} : {len(cl_residues_probs[key].keys())} residues{os.linesep}"
            cl_string += f"rank\tresid\tprobability{os.linesep}"
            sorted_probs = sorted(
                cl_residues_probs[key].items(), key=lambda x: x[1], reverse=True
            )
            for pair_idx, pair in enumerate(sorted_probs, start=1):
                cl_string += f"{pair_idx}\t{pair[0]}\t{pair[1]:.3f}{os.linesep}"
            cl_string += os.linesep
            wfile.write(cl_string)


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
    cl_dict = get_clustering_dict(clusters, entries)
    write_clusters(cl_dict, cl_filename)
    # write clustered residues
    res_filename = "clustered_residues.out"
    cl_residues, cl_residues_probs = get_residue_dict(cl_dict, interface_dict)
    write_residues(cl_residues, res_filename)
    res_probs_filename = "clustered_residues_probs.out"
    write_residues_probs(cl_residues_probs, res_probs_filename)
    # write time
    elap_time = round((time.time() - start_time), 3)
    log.info(f"Clustering performed in {elap_time} seconds")
    return cl_residues, cl_residues_probs
