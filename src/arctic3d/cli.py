"""Main CLI."""
import argparse
import sys
import time
from pathlib import Path

from arctic3d import log
from arctic3d.modules.blast import run_blast
from arctic3d.modules.clustering import filter_clusters
from arctic3d.modules.cluster_interfaces import cluster_interfaces
from arctic3d.modules.input import Input
from arctic3d.modules.interface import (
    get_interface_residues,
    read_interface_residues,
)
from arctic3d.modules.output import (
    get_init_message,
    make_output,
    create_output_folder,
    setup_output_folder,
)
from arctic3d.modules.pdb import get_best_pdb
from arctic3d.modules.sequence import to_fasta
from arctic3d.modules.log import add_log_for_CLI


argument_parser = argparse.ArgumentParser()
argument_parser.add_argument(
    "input_arg",
    help="",
)

argument_parser.add_argument(
    "--db",
    help="",
)

argument_parser.add_argument(
    "--interface_file",
    help="input interface file",
)

argument_parser.add_argument(
    "--out_partner",
    help="set of comma-separated partner IDs to exclude from the search",
)

argument_parser.add_argument(
    "--out_pdb",
    help="set of comma-separated pdb IDs to exclude from the search",
)

argument_parser.add_argument(
    "--pdb_to_use",
    help="pdb file to be used",
)

argument_parser.add_argument(
    "--chain_to_use",
    help="the chain to be used",
)

argument_parser.add_argument(
    "--run_dir",
    help="directory where to store the run",
)

argument_parser.add_argument(
    "--interface_data",
    help=".json file containing the interface data",
)

argument_parser.add_argument(
    "--pdb_data",
    help=".json file containing the pdb data",
)

argument_parser.add_argument(
    "--full",
    help="consider full uniprot-pdb-chain information in the retrieval",
    action="store_true",
)

argument_parser.add_argument(
    "--ligand",
    help="retrieve ligand binding residues",
    default="no",
    choices=["yes", "no", "both"],
)

argument_parser.add_argument(
    "--threshold",
    help="Threshold for clustering",
    type=float,
    required=False,
    default=0.866,
)

argument_parser.add_argument(
    "--linkage_strategy",
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
    "--min_clust_size",
    help="Minimum number of residues in clusters",
    type=int,
    required=False,
    default=0,
)

argument_parser.add_argument(
    "--int_cov_cutoff",
    help="Interface coverage cutoff (%%)",
    type=float,
    required=False,
    default=0.7,
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
    return main_func(**vars(cmd))


def maincli():
    """Execute main client."""
    return cli(argument_parser, main)


def main(
    input_arg,
    db,
    interface_file,
    out_partner,
    out_pdb,
    pdb_to_use,
    chain_to_use,
    run_dir,
    interface_data,
    pdb_data,
    full,
    ligand,
    linkage_strategy,
    threshold,
    min_clust_size,
    int_cov_cutoff,
    log_level="DEBUG",
):
    """Main function."""
    init_message = get_init_message()
    log.info(init_message)
    st_time = time.time()
    inp = Input(input_arg)
    input_files = {}
    # retrieve uniprot information
    if inp.is_fasta():
        input_files["fasta"] = Path(inp.arg)
        uniprot_id = run_blast(input_files["fasta"], db)
    if inp.is_uniprot():
        uniprot_id = inp.arg.upper()
    if inp.is_pdb():
        input_files["pdb"] = Path(inp.arg)
        if not interface_file:
            fasta_f = to_fasta(input_files["pdb"], temp=False)
            uniprot_id = run_blast(fasta_f.name, db)
        else:
            input_files["interface_file"] = Path(interface_file)
            uniprot_id = None

    # create output folder
    run_dir_path = create_output_folder(run_dir, uniprot_id)
    # configure logging
    log_file = Path(run_dir_path, "arctic3d.log")
    add_log_for_CLI(log, log_level, log_file)

    log.info(f"Target UNIPROTID: {uniprot_id}")

    # save json files
    if interface_data:
        input_files["interface_data"] = Path(interface_data)
    if pdb_data:
        input_files["pdb_data"] = Path(pdb_data)

    input_files = setup_output_folder(
        run_dir=run_dir_path, input_files=input_files
    )

    # retrieve interfaces.
    if "interface_file" in input_files:
        log.info(f"input interface file {interface_file}")
        interface_residues = read_interface_residues(
            input_files["interface_file"]
        )
    else:
        if interface_file:
            log.warning(
                "input interface file submitted without pdb. It will be"
                " ignored."
            )
        if interface_data:
            int_data_path = input_files["interface_data"]
        else:
            int_data_path = None
        interface_residues = get_interface_residues(
            uniprot_id=uniprot_id,
            out_partner_string=out_partner,
            out_pdb_string=out_pdb,
            full=full,
            ligand=ligand,
            interface_data=int_data_path,
        )

    log.info(f"Interface Residues: {interface_residues}")

    if interface_residues:
        # retrieve pdb file
        if inp.is_pdb():
            # interfaces will be filtered later
            pdb_f, filtered_interfaces = input_files["pdb"], None
            if not interface_file:
                log.warning(
                    (
                        "Input pdb file submitted without interface file. "
                        "This assumes the pdb is coherent with the "
                        "corresponding uniprot numbering."
                    )
                )
        else:
            if pdb_data:
                pdb_data_path = input_files["pdb_data"]
            else:
                pdb_data_path = None
            # get best pdb
            pdb_f, cif_f, filtered_interfaces = get_best_pdb(
                uniprot_id=uniprot_id,
                interface_residues=interface_residues,
                pdb_to_use=pdb_to_use,
                chain_to_use=chain_to_use,
                pdb_data=pdb_data_path,
                int_cov_cutoff=int_cov_cutoff,
            )

        if pdb_f is None:
            log.error(
                "Could not retrieve a valid PDB for the target, please provide"
                " one as the main input argument."
            )
            sys.exit()

        log.info(f"PDB file: {pdb_f}")

        # cluster interfaces
        if filtered_interfaces:
            interface_dict = filtered_interfaces
        else:
            interface_dict = interface_residues
        cl_dict, cl_residues, cl_residues_probs = cluster_interfaces(
            interface_dict=interface_dict,
            pdb_path=pdb_f,
            linkage_strategy=linkage_strategy,
            threshold=threshold,
            int_cov_cutoff=int_cov_cutoff,
        )

        log.info(f"Clustered interfaces {cl_dict}")
        log.info(f"Clustered interface residues: {cl_residues}")
        if min_clust_size > 0:
            log.info(
                f"Excluding clusters with less than {min_clust_size} residues"
            )
            cl_dict, cl_residues, cl_residues_probs = filter_clusters(
                cl_dict, cl_residues, cl_residues_probs, min_clust_size
            )

        make_output(
            interface_residues=interface_residues,
            pdb_f=pdb_f,
            cl_dict=cl_dict,
            cl_residues=cl_residues,
            cl_residues_probs=cl_residues_probs,
        )
    else:
        log.info("No interfaces found.")

    log.info(
        f"arctic3d run completed in {(time.time() - st_time):.2f} seconds."
    )

    # check if there's at least one interface
    if Path("clustered_interfaces.out").is_file() is False:
        return 255

    return 0


if __name__ == "__main__":
    sys.exit(maincli())
