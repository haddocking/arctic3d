"""Main CLI."""
import argparse
import sys

from arctic3d.modules.blast import blast_seq
from arctic3d.modules.interface import get_interface_residues
from arctic3d.modules.output import make_output
from arctic3d.modules.pdb import download_pdb, get_best_pdb
from arctic3d.modules.sequence import load_seq

argument_parser = argparse.ArgumentParser()
argument_parser.add_argument(
    "fasta_file",
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


def main(fasta_file):
    """Main function."""

    fasta_seq = load_seq(fasta_file)

    uniprot_id = blast_seq(fasta_seq)

    best_pdb = get_best_pdb(uniprot_id)

    _ = download_pdb(best_pdb)

    interface_residues = get_interface_residues(uniprot_id)

    _ = make_output(best_pdb, interface_residues)


if __name__ == "__main__":
    sys.exit(maincli())
