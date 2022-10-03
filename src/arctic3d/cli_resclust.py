"""
Residue-based clustering.

Given a pdb and a set of residues, you can create different clusters of
residues.

USAGE::

    arctic3d_resclust ./example/1ppe_E.pdb --residue_list 29,30,31,49,50,51 --threshold=20.0 --chain=E

Input arguments:
`residue_list` : the comma-separated list of residue IDs
`threshold` : the number to be used as threshold for hierarchical clustering
`chain` : the chain ID to be used
"""
import argparse
import logging
import sys

import MDAnalysis as mda
from scipy.spatial.distance import pdist

from arctic3d.modules.clustering import cluster_similarity_matrix, get_clustering_dict
from arctic3d.modules.input import Input

log = logging.getLogger("arctic3dlog")
ch = logging.StreamHandler()
formatter = logging.Formatter(
    " [%(asctime)s %(module)s:L%(lineno)d %(levelname)s] %(message)s"
)
ch.setFormatter(formatter)
log.addHandler(ch)

THRESHOLD = 10.0

argument_parser = argparse.ArgumentParser()
argument_parser.add_argument(
    "input_arg",
    help="Input PDB file",
)

argument_parser.add_argument(
    "--residue_list", help="List of residues to cluster", required=True
)

argument_parser.add_argument(
    "--threshold",
    help="Threshold (in angstroms) for clustering",
    type=float,
    required=False,
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


def main(input_arg, residue_list, threshold, chain):
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

    log.info(f"resids_list {resids_list}")
    str_resids_list = [str(res) for res in resids_list]
    sel_residues = f"name CA {chain_str}and resid {' '.join(str_resids_list)}"

    u = mdu.select_atoms(sel_residues)
    log.info(f"retrieved residues: {u.resids}")

    n_chains = u.n_segments
    if n_chains != 1:
        log.error(f"Number of consistent segments ({n_chains}) != 1. Aborting.")
        sys.exit(1)

    # do the clustering
    cutoff = THRESHOLD
    if threshold:
        cutoff = threshold
    distmap = pdist(u.positions)
    clusters = cluster_similarity_matrix(distmap, resids_list, threshold=cutoff)
    cl_dict = get_clustering_dict(clusters, resids_list)
    for el in cl_dict.keys():
        log.info(
            f"cluster {el}, residues {' '.join([str(res) for res in cl_dict[el]])}"
        )


if __name__ == "__main__":
    sys.exit(maincli())
