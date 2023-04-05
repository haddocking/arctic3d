"""
Residue-based clustering.

Given a pdb and a set of residues, you can create different clusters of
residues.

USAGE::

    arctic3d_resclust ./example/1ppe_E.pdb \
        --residue_list 29,30,31,49,50,51 \
        --threshold=20.0 \
        --chain=E

Input arguments:
    `residue_list` : the comma-separated list of residue IDs.

    `threshold` : the number to be used as threshold for hierarchical
        clustering. If `criterion` is `maxclust`, this is the maximum number
        of clusters.

    `chain` : the chain ID to be used.

    `linkage` : the linkage strategy.

    `criterion` : the criterion to extract the clusters.
"""
import argparse
import sys

import MDAnalysis as mda
from scipy.spatial.distance import pdist

from arctic3d import log
from arctic3d.modules.clustering import (
    cluster_similarity_matrix,
    get_clustering_dict,
)
from arctic3d.modules.input import Input


argument_parser = argparse.ArgumentParser()
argument_parser.add_argument(
    "input_arg",
    help="Input PDB file",
)

argument_parser.add_argument(
    "--residue_list",
    help="List of (comma-separated) residues to cluster",
    required=True,
)

argument_parser.add_argument(
    "--threshold",
    help="Threshold (in angstroms) for clustering",
    type=float,
    required=False,
    default=15.0,
)

argument_parser.add_argument(
    "--criterion",
    help="Criterion for clustering",
    type=str,
    required=False,
    choices=["distance", "maxclust"],
    default="distance",
)

argument_parser.add_argument(
    "--linkage",
    help="Linkage strategy for clustering",
    type=str,
    required=False,
    choices=[
        "average",
        "single",
        "complete",
        "median",
        "centroid",
        "ward",
        "weighted",
    ],
    default="average",
)

argument_parser.add_argument(
    "--chain", help="Segment ID to be considered", required=False
)


def load_args(arguments):
    """
    Load argument parser.

    Parameters
    ----------
    arguments : argparse.ArgumentParser
        Argument parser.

    Returns
    -------
    cmd : argparse.Namespace
        Parsed command-line arguments.

    """
    return arguments.parse_args()


def cli(arguments, main_func):
    """
    Command-line interface entry point.

    Parameters
    ----------
    arguments : argparse.ArgumentParser
        Argument parser.
    main_func : function
        Main function.

    """
    cmd = load_args(arguments)
    main_func(**vars(cmd))


def maincli():
    """Execute main client."""
    cli(argument_parser, main)


def main(input_arg, residue_list, chain, threshold, linkage, criterion):
    """Main function."""
    log.setLevel("INFO")

    # check input
    inp = Input(input_arg)
    if not inp.is_pdb:
        log.error("Input must be a pdb file")
        sys.exit(1)

    # read pdb
    try:
        mdu = mda.Universe(inp.arg)
    except ValueError:
        log.error(f"Unable to read input PDB file {inp.arg}")
        sys.exit(1)

    # extract atoms
    if chain:
        chain_str = f"and chainID {chain} "
    else:
        chain_str = ""

    try:
        resids_list = [int(el) for el in residue_list.strip().split(",")]
    except ValueError:
        log.error(f"Malformed input residue_list {residue_list}")
        sys.exit(1)

    log.info(f"input residue_list {resids_list}")
    str_resids_list = [str(res) for res in resids_list]
    sel_residues = f"name CA {chain_str}and resid {' '.join(str_resids_list)}"

    u = mdu.select_atoms(sel_residues)
    unique_sorted_resids = u.resids
    log.info(f"retrieved residues: {unique_sorted_resids}")

    n_chains = u.n_segments
    if n_chains != 1:
        log.error(f"Number of consistent segments ({n_chains}) != 1.Aborting.")
        sys.exit(1)

    # do the clustering
    if criterion == "maxclust":
        threshold = int(threshold)
    log.info(
        f"Clustering distance matrix with linkage {linkage}, threshold"
        f" {threshold}, and criterion {criterion}"
    )
    distmap = pdist(u.positions)
    clusters = cluster_similarity_matrix(
        distmap,
        unique_sorted_resids,
        threshold=threshold,
        linkage_strategy=linkage,
        crit=criterion,
    )

    cl_dict = get_clustering_dict(clusters, unique_sorted_resids)
    for el in cl_dict.keys():
        log.info(
            f"cluster {el}, residues"
            f" {' '.join([str(res) for res in cl_dict[el]])}"
        )


if __name__ == "__main__":
    sys.exit(maincli())
