"""Main CLI."""
import argparse
import logging
import sys
from pathlib import Path

from arctic3d.modules.blast import run_blast
from arctic3d.modules.sequence import to_fasta

log = logging.getLogger("arctic3dlog")
ch = logging.StreamHandler()
formatter = logging.Formatter(
    " [%(asctime)s %(module)s:L%(lineno)d %(levelname)s] %(message)s"
)
ch.setFormatter(formatter)
log.addHandler(ch)


argument_parser = argparse.ArgumentParser()
argument_parser.add_argument("input_pdbs", nargs="+", default=[])

argument_parser.add_argument(
    "--db",
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


def main(input_pdbs, db):
    """Main function."""
    log.setLevel("DEBUG")

    with open("uniprot.csv", "w") as fh:
        for pdb in input_pdbs:
            fasta_f = to_fasta(pdb, temp=True)
            uniprot_id = run_blast(fasta_f.name, db)
            Path(fasta_f.name).unlink()
            log.info(f"pdb: {pdb} uniprot_id: {uniprot_id}")
            fh.write(f"{pdb},{uniprot_id}\n")


if __name__ == "__main__":
    sys.exit(maincli())
