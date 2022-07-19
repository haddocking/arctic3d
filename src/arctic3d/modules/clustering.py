"""Clustering module."""

import logging
import numpy as np
import pandas as pd
import os
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram
import matplotlib.pyplot as plt


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
    print(f"nligands {nligands} vs int_nligands {int_nligands}")
    if abs(nligands - int_nligands) > 0.00001:
        raise Exception(f"npairs {int_matrix.shape[0]}: interface matrix should be a 1D condensed distance matrix")
    ligand_names = [int_matrix.iloc[0,0]]
    for lig in np.unique(int_matrix.iloc[:,1]):
        ligand_names.append(lig)
    print(f"ordered ligand names {ligand_names}")
    return int_matrix.iloc[:,2], ligand_names

def cluster_distance_matrix(int_matrix, entries, plot=None):
    """
    Does the clustering.
    
    Parameters
    ----------
    
    Returns
    -------

    """
    print(f"creating dendrogram")
    Z = linkage(int_matrix, LINKAGE)
    if plot:
        dendrogram_figure_filename = plot + "_" + LINKAGE + ".png"
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
    print(f"clusters = {clusters}")
    return clusters

def write_clusters(clusters, ligands, cl_filename):
    """writes clusters to file."""
    print(f"Writing clusters to file {cl_filename}")
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


def interface_clustering(matrix_filename):
    """clusters the interface matrix."""
    # read matrix
    if os.path.exists(matrix_filename):
        int_matrix, entries = read_int_matrix(matrix_filename)
    else:
        raise Exception(f"input path {matrix_filename} does not exist!")
    # cluster matrix
    clusters = cluster_distance_matrix(int_matrix, entries)
    # write clusters
    cl_filename = "interface-clustering_" + LINKAGE + "_thr-" + str(THRESHOLD) + ".out"
    write_clusters(clusters, entries, cl_filename)

interface_clustering("./interface_matrix.txt")