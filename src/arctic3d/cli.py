"""Main CLI."""
import argparse
import logging
import sys
from pathlib import Path

from arctic3d.modules.blast import run_blast
from arctic3d.modules.cluster_interfaces import cluster_interfaces

# from arctic3d.modules.geometry import cluster_interface
from arctic3d.modules.input import Input
from arctic3d.modules.interface import get_interface_residues, read_interface_residues
from arctic3d.modules.output import make_output, setup_output_folder
from arctic3d.modules.pdb import get_best_pdb
from arctic3d.modules.sequence import to_fasta

# from arctic3d.modules.sequence import load_seq

log = logging.getLogger("arctic3dlog")
ch = logging.StreamHandler()
formatter = logging.Formatter(
    " [%(asctime)s %(module)s:L%(lineno)d %(levelname)s] %(message)s"
)
ch.setFormatter(formatter)
log.addHandler(ch)


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
    help="",
)

argument_parser.add_argument(
    "--out_uniprot",
    help="",
)

argument_parser.add_argument(
    "--out_pdb",
    help="",
)

argument_parser.add_argument(
    "--pdb_to_use",
    help="",
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
    default=False,
)

argument_parser.add_argument(
    "--ligand",
    help="retrieve ligand binding residues",
    default="no",
    choices=["yes", "no", "both"]
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


def main(
    input_arg,
    db,
    interface_file,
    out_uniprot,
    out_pdb,
    pdb_to_use,
    chain_to_use,
    run_dir,
    interface_data,
    pdb_data,
    full,
    ligand,
):
    """Main function."""
    log.setLevel("DEBUG")

    inp = Input(input_arg)
    input_files = {}
    # retrieve uniprot information
    if inp.is_fasta():
        input_files["fasta"] = Path(inp.arg)
        uniprot_id = run_blast(input_files["fasta"], db)
    if inp.is_uniprot():
        uniprot_id = inp.arg
    if inp.is_pdb():
        input_files["pdb"] = Path(inp.arg)
        if not interface_file:
            fasta_f = to_fasta(input_files["pdb"], temp=False)
            uniprot_id = run_blast(fasta_f.name, db)
        else:
            input_files["interface_file"] = Path(interface_file)
            uniprot_id = None

    # save json files
    if interface_data:
        input_files["interface_data"] = Path(interface_data)
    if pdb_data:
        input_files["pdb_data"] = Path(pdb_data)

    log.info(f"Target UNIPROTID: {uniprot_id}")

    input_files = setup_output_folder(uniprot_id, input_files, run_dir)

    # retrieve interfaces.
    if interface_file:
        log.info(f"input interface file {interface_file}")
        interface_residues = read_interface_residues(input_files["interface_file"])
    else:
        if interface_data:
            interface_residues = get_interface_residues(
                uniprot_id, out_uniprot, out_pdb, input_files["interface_data"], full
            )
        else:
            interface_residues = get_interface_residues(
                uniprot_id, out_uniprot, out_pdb, full, ligand
            )

    log.info(f"Interface Residues: {interface_residues}")

    if interface_residues:
        # retrieve pdb file
        if inp.is_pdb():
            # interfaces will be filtered later
            pdb_f, filtered_interfaces = input_files["pdb"], None
            if not interface_file:
                log.warning(
                    """Input pdb file submitted without interface file. This assumes the pdb is coherent with the corresponding uniprot numbering."""
                )
        else:
            if pdb_data:
                pdb_f, filtered_interfaces = get_best_pdb(
                    uniprot_id,
                    interface_residues,
                    pdb_to_use,
                    chain_to_use,
                    input_files["pdb_data"],
                )
            else:
                pdb_f, filtered_interfaces = get_best_pdb(
                    uniprot_id, interface_residues, pdb_to_use, chain_to_use
                )

        if pdb_f is None:
            log.error(
                "Could not retrieve a valid PDB for the target, please provide one using the --pdb option"
            )
            sys.exit()

        log.info(f"PDB file: {pdb_f}")

        # cluster interfaces
        if filtered_interfaces:
            cl_ints, cl_residues, cl_residues_probs = cluster_interfaces(
                filtered_interfaces, pdb_f
            )
        else:
            cl_ints, cl_residues, cl_residues_probs = cluster_interfaces(
                interface_residues, pdb_f
            )
        log.info(f"Clustered interfaces {cl_ints}")
        log.info(f"Clustered interface residues: {cl_residues}")

        make_output(pdb_f, cl_ints, cl_residues, cl_residues_probs)
    else:
        log.info("No interfaces found.")


if __name__ == "__main__":
    sys.exit(maincli())
