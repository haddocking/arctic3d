"""
Get HADDOCK restraints from arctic3d data.

Given two arctic3d run folders, this script will generate a HADDOCK
restraints file for the interface residues present in the
clustered_residues_probs.out file.

USAGE::

    arctic3d_restraints --r1 ./arctic3d_run1 --r2 ./arctic3d_run2

Use the ch1 and ch2 parameters if you want to specify different chains to be
used for the restraints (default is A for r1 and B for r2)::

    arctic3d_restraints --r1 ./arctic3d_run1 --r2 ./arctic3d_run2 \
        --ch1=X --ch2=Y

Remember that the chain IDs must be consistent with the chain present in your
PDB file.

Use the run_dir parameter if you want to specify a specific output directory::

    arctic3d_restraints --r1 ./arctic3d_run1 --r2 ./arctic3d_run2 \
        --run_dir=arctic3d-restraints-example

Use the prob_threshold parameter to specify the probability threshold for
the interface residues::

    arctic3d_restraints --r1 ./arctic3d_run1 --r2 ./arctic3d_run2 \
        --prob_threshold=0.5

This will consider only residues with a probability of being in the interface
higher than 0.5 (for each cluster).
"""
import argparse
import logging
import os
import shutil
import sys
import time
from pathlib import Path
import tarfile

from arctic3d.modules.output import read_residues_probs, setup_output_folder

LOGNAME = f"arctic3d_restraints_{os.getpid()}.log"
logging.basicConfig(filename=LOGNAME)
log = logging.getLogger(LOGNAME)
ch = logging.StreamHandler()
formatter = logging.Formatter(
    " [%(asctime)s %(module)s:L%(lineno)d %(levelname)s] %(message)s"
)
ch.setFormatter(formatter)
log.addHandler(ch)

argument_parser = argparse.ArgumentParser()
argument_parser.add_argument(
    "--r1",
    required=True,
    help="Directory of run1",
)

argument_parser.add_argument(
    "--r2",
    required=True,
    help="Directory of run2",
)

argument_parser.add_argument(
    "--ch1",
    default="A",
    help="chain ID for run1 residues",
)

argument_parser.add_argument(
    "--ch2",
    default="B",
    help="chain ID for run2 residues",
)

argument_parser.add_argument(
    "--run_dir",
    help="directory where to store the run",
    default="arctic3d-restraints",
)

argument_parser.add_argument(
    "--prob_threshold",
    help="probability threshold for interface residues",
    type=float,
    default=0.3,
)


def filter_residues_probs(residues_probs, prob_threshold):
    """Filter residues_probs based on the probability threshold.

    Parameters
    ----------
    residues_probs : dict
        Dictionary with the residues probabilities.
    prob_threshold : float
        Probability threshold.

    Returns
    -------
    filtered_residues_probs : dict of lists
        Dictionary with the clustered residues filtered by probability.
    """
    filtered_residues_probs = {}
    for cluster, residues in residues_probs.items():
        res_list = []
        for residue, prob in residues.items():
            if prob >= prob_threshold:
                res_list.append(residue)
        # appending sorted list (to have a well formatted output)
        filtered_residues_probs[cluster] = sorted(res_list)
    return filtered_residues_probs


def assign_line(resid, chain):
    """return assign line."""
    ass_line = "assign ( resid "
    ass_line += f"{resid} and segid {chain})"
    ass_line += f"{os.linesep}       ({os.linesep}"
    return ass_line


def generate_restraints(residues1, residues2, ch1, ch2, ambig_fname):
    """
    Generate act-act.sh restraint file.

    Parameters
    ----------
    residues1 : list
        List of residues for the first partner.
    residues2 : list
        List of residues for the second partner.
    ch1 : str
        Chain ID for the first partner.
    ch2 : str
        Chain ID for the second partner.
    ambig_fname : str
        Name of the output file.
    """
    with open(ambig_fname, "w") as ambig_file:
        ambig_file.write(
            f"!HADDOCK AIR restraints for 1st partner{os.linesep}!{os.linesep}"
        )
        for res in residues1:
            ambig_file.write(assign_line(res, ch1))
            res2_group = [
                f"        ( resid {res2} and segid {ch2})"
                for res2 in residues2
            ]
            res2_string = f"{os.linesep}     or{os.linesep}".join(res2_group)
            ambig_file.write(f"{res2_string}{os.linesep}")
            ambig_file.write(f"       )  2.0 2.0 0.0{os.linesep}!{os.linesep}")
        # second partner
        ambig_file.write(
            f"!{os.linesep}!HADDOCK AIR restraints for "
            f"2nd partner{os.linesep}!{os.linesep}"
        )
        for res in residues2:
            ambig_file.write(assign_line(res, ch2))
            res1_group = [
                f"        ( resid {res1} and segid {ch1})"
                for res1 in residues1
            ]
            res1_string = f"{os.linesep}     or{os.linesep}".join(res1_group)
            ambig_file.write(f"{res1_string}{os.linesep}")
            ambig_file.write(f"       )  2.0 2.0 0.0{os.linesep}!{os.linesep}")
    return


def compress_tbl_files(ambig_fnames, out_tgz):
    """
    Compress restraints in a tbl.tgz file.

    Parameters
    ----------
    ambig_fnames : list of Path
        List of restraint files.
    out_tgz : str
        Name of the output tgz file.
    """
    # compress restraints in a tbl.tgz file
    log.info(f"Compressing restraints into {out_tgz}")
    tgz = tarfile.open(out_tgz, "w:gz")
    for name in ambig_fnames:
        tgz.add(name)
    tgz.close()

    return


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


def main(r1, r2, ch1, ch2, run_dir, prob_threshold=0.5):
    """Main function."""
    log.setLevel("INFO")
    start_time = time.time()
    log.info("Starting arctic3d_restraints")

    # checking if r1 and r2 exists
    if not os.path.exists(r1):
        log.error(f"Could not find {r1}")
        sys.exit(1)
    if not os.path.exists(r2):
        log.error(f"Could not find {r2}")
        sys.exit(1)

    # checking if r1 and r2 contain the clustered_residues_probs.out file
    r1_res_fname = Path(r1, "clustered_residues_probs.out")
    r2_res_fname = Path(r2, "clustered_residues_probs.out")
    if not os.path.exists(r1_res_fname):
        log.error(f"Could not find clustered_residues_probs.out in {r1}")
        sys.exit(1)
    if not os.path.exists(r2_res_fname):
        log.error(f"Could not find clustered_residues_probs.out in {r2}")
        sys.exit(1)

    # Setting up output folder
    input_files = {"r1_res_fname": r1_res_fname, "r2_res_fname": r2_res_fname}
    log.info(f"Input files are {input_files}")
    input_files = setup_output_folder(None, input_files, run_dir)

    # read and filter probabilities
    r1_residues_probs = read_residues_probs(input_files["r1_res_fname"])
    r2_residues_probs = read_residues_probs(input_files["r2_res_fname"])
    r1_residues = filter_residues_probs(r1_residues_probs, prob_threshold)
    r2_residues = filter_residues_probs(r2_residues_probs, prob_threshold)
    log.info(f"{r1} residues = {r1_residues}")
    log.info(f"{r2} residues = {r2_residues}")

    # creating restraints from residues
    ambig_fnames = []
    tot_nrestrs = len(r1_residues) * len(r2_residues)
    n_ambig = 0
    log.info(f"Creating {tot_nrestrs} restraints")
    for cl1, residues1 in r1_residues.items():
        for cl2, residues2 in r2_residues.items():
            ambig_fname = f"ambig_{n_ambig}.tbl"
            ambig_fnames.append(ambig_fname)
            log.info(
                f"Creating {ambig_fname} restraint file by"
                "coupling {cl1} (r1) and {cl2} (r2)"
            )
            generate_restraints(residues1, residues2, ch1, ch2, ambig_fname)
            n_ambig += 1

    compress_tbl_files(ambig_fnames, out_tgz="ambig.tbl.tgz")

    elap_time = round((time.time() - start_time), 3)
    log.info(f"arctic3d_restraints run took {elap_time} seconds")

    # copying log file to the run folder (if possible)
    try:
        shutil.move(f"../{LOGNAME}", LOGNAME)
    except FileNotFoundError as e:
        log.warning(f"Could not find log file: {e}")
