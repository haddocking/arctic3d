"""Clustering module."""

import logging
import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage

from arctic3d.modules.interface_matrix import read_int_matrix

log = logging.getLogger("arctic3d.log")


def plot_dendrogram(
    linkage_matrix, entries, filename, threshold, max_entries=50
):
    """
    Plots the dendrogram.

    Parameters
    ----------
    linkage matrix : np.array (4 x nentries)
        linkage matrix
    entries : list
        list of interface names
    filename : str or Path
        plot filename
    threshold : float
        threshold for coloring
    max_entries : int
        maximum number of entries to plot
    """
    plt.figure(dpi=200)
    truncate_mode = None
    p = len(entries)
    if len(entries) > max_entries:
        log.info("High number of entries, truncating dendrogram...")
        truncate_mode = "lastp"
        p = max_entries
    dendrogram(
        linkage_matrix,
        color_threshold=threshold,
        labels=entries,
        truncate_mode=truncate_mode,
        p=p,
        orientation="right",
    )
    plt.ylabel("Interface Names")
    plt.xlabel("Dissimilarity")
    plt.title("ARCTIC3D dendrogram")
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()


def cluster_similarity_matrix(
    int_matrix,
    entries,
    linkage_strategy="average",
    threshold=0.866,
    crit="distance",
    plot=False,
):
    """
    Does the clustering.

    Parameters
    ----------
    int_matrix : np.array
        1D condensed interface similarity matrix
    entries : list
        names of the ligands
    linkage_strategy : str
        linkage strategy for clustering
    threshold : float
        threshold for clustering
    plot : bool
        if True, plot the dendrogram

    Returns
    -------
    clusters : list
        list of clusters ID, each one associated to an entry
    """
    log.info(
        f"Clustering with linkage {linkage_strategy} and threshold {threshold}"
    )
    Z = linkage(int_matrix, linkage_strategy)
    if plot:
        dendrogram_figure_filename = "dendrogram_" + linkage_strategy + ".png"
        plot_dendrogram(Z, entries, dendrogram_figure_filename, threshold)
    # clustering
    clusters = fcluster(Z, t=threshold, criterion=crit)
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
        *example* {
            1 : [
            'interface_1', 'interface_3'
                ] ,
            2 : [
            'interface_2'
            ],
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
        *example* { 1 : [1,2,3,5,6,8] , 2 : [29,30,31], ... }
    cl_residues_probs : dict of dicts
        dictionary of probabilities for clustered residues
        *example* { 1 : {1:0.7, 2:0.2, 3:0.4 ...}, ... }
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


def interface_clustering(
    interface_dict, matrix_filename, linkage_strategy, threshold
):
    """
    Clusters the interface matrix.

    Parameters
    ----------
    interface_dict : dict
        dictionary of all the interfaces (each one with its uniprot ID as key)
    matrix_filename : str or Path
        input interface matrix
    linkage_strategy : str
        linkage strategy for clustering
    threshold : float
        threshold for clustering

    Returns
    -------
    cl_dict : dict
        dictionary of clustered interfaces
    cl_residues : dict
        dictionary of clustered residues
    cl_residues_probs : dict of dicts
        dictionary of probabilities for clustered residues
    """
    start_time = time.time()
    log.info("Clustering interface matrix")

    # check if there's only a single interface
    if len(interface_dict) == 1:
        clusters = [1]
        entries = list(interface_dict.keys())  # the only entry
    else:
        int_matrix, entries = read_int_matrix(matrix_filename)  # read matrix
        # cluster matrix.
        clusters = cluster_similarity_matrix(
            int_matrix, entries, linkage_strategy, threshold, plot=True
        )

    # get clustering dictionary and clustered_residues
    cl_dict = get_clustering_dict(clusters, entries)
    cl_residues, cl_residues_probs = get_residue_dict(cl_dict, interface_dict)

    # write time
    elap_time = round((time.time() - start_time), 3)
    log.info(f"Clustering performed in {elap_time} seconds")
    log.info(f"Clustering produced {len(cl_dict)} clusters")
    return cl_dict, cl_residues, cl_residues_probs
