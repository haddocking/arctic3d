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
from arctic3d.modules.output import setup_output_folder
from arctic3d.modules.pdb import get_best_pdb, output_pdb
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


def main(input_arg, db, interface_file, out_uniprot, out_pdb, pdb_to_use):
    """Main function."""
    log.setLevel("DEBUG")

    inp = Input(input_arg)
    input_files = []
    # retrieve uniprot information
    if inp.is_fasta():
        fasta_f = Path(inp.arg)
        uniprot_id = run_blast(fasta_f, db)
        input_files.append(fasta_f)
    if inp.is_uniprot():
        uniprot_id = inp.arg
    if inp.is_pdb():
        pdb_f = Path(inp.arg)
        input_files.append(pdb_f)
        if not interface_file:
            fasta_f = to_fasta(pdb_f, temp=False)
            uniprot_id = run_blast(fasta_f.name, db)
        else:
            interface_f = Path(interface_file)
            uniprot_id = None
            input_files.append(Path(interface_file))

    log.info(f"Target UNIPROTID: {uniprot_id}")

    setup_output_folder(uniprot_id, input_files)

    # retrieve interfaces
    if interface_file:
        log.info(f"input interface file {interface_file}")
        interface_residues = read_interface_residues(
            Path("input_data", interface_f.name)
        )
    else:
        interface_residues = get_interface_residues(uniprot_id, out_uniprot, out_pdb)

    log.info(f"Interface Residues: {interface_residues}")

    if interface_residues:
        # retrieve pdb file
        if inp.is_pdb():
            # interfaces will be filtered later
            pdb_f, filtered_interfaces = Path("input_data", pdb_f.name), None
            if not interface_file:
                log.warning(
                    """Input pdb file submitted without interface file. This assumes the pdb is coherent with the corresponding uniprot numbering."""
                )
        else:
            print(f"pdb_to_use {pdb_to_use}")
            pdb_f, filtered_interfaces = get_best_pdb(
                uniprot_id, interface_residues, pdb_to_use
            )

        if pdb_f is None:
            log.error(
                "Could not retrieve a valid PDB for the target, please provide one using the --pdb option"
            )
            sys.exit()

        log.info(f"PDB file: {pdb_f}")

        # cluster interfaces
        if filtered_interfaces:
            clustered_interface_residues, cl_residues_probs = cluster_interfaces(
                filtered_interfaces, pdb_f
            )
        else:
            clustered_interface_residues, cl_residues_probs = cluster_interfaces(
                interface_residues, pdb_f
            )

        log.info(f"Clustered interface residues: {clustered_interface_residues}")

        output_pdb(pdb_f, cl_residues_probs)
    else:
        log.info("No interfaces found.")


if __name__ == "__main__":
    sys.exit(maincli())
