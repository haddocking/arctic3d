import argparse
import logging
import os
import shlex
import subprocess
import time
from pathlib import Path

import pandas as pd

logging.basicConfig(
    format=" %(funcName)s:L%(lineno)d %(levelname)s - %(message)s",
    level=logging.INFO,
    datefmt="%m/%d/%Y %I:%M:%S %p",
)

INT_FILENAME = "clustered_interfaces.out"


def get_arctic_io(filename):
    """creates a dictionary of arctic3d input-output."""
    bm5_uniprot = pd.read_csv(filename)
    arctic_io = {}
    for n in range(bm5_uniprot.shape[0]):
        complex_pdb = bm5_uniprot["complex"].iloc[n]
        receptor_pdb = bm5_uniprot["receptor"].iloc[n]
        receptor_uniprot = bm5_uniprot["uniprot_receptor"].iloc[n]
        ligand_pdb = bm5_uniprot["ligand"].iloc[n]
        ligand_uniprot = bm5_uniprot["uniprot_ligand"].iloc[n]

        arctic_io[f"{complex_pdb}-1"] = {
            "complex_pdb": complex_pdb,
            "receptor_pdb": receptor_pdb,
            "self_uniprot_id": receptor_uniprot,
            "paired_uniprot_id": ligand_uniprot,
            "arctic_output": {},
        }

        arctic_io[f"{complex_pdb}-2"] = {
            "complex_pdb": complex_pdb,
            "receptor_pdb": ligand_pdb,
            "self_uniprot_id": ligand_uniprot,
            "paired_uniprot_id": receptor_uniprot,
            "arctic_output": {},
        }

    return arctic_io


def run_arctic(cmd):
    """runs arctic3d command and returns the status"""
    p = subprocess.run(
        shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    # check if failed
    if p.returncode != 0:
        logging.info(f"warning: ARCTIC failed with cmd {cmd}")
        response = "FAILED"
    else:
        response = "SUCCESS"
    return response


def read_interface_file(interface_filename):
    """Reads interface file"""
    interfaces = {}
    with open(interface_filename, "r") as ofile:
        for ln in ofile:
            splt_ln = ln.split()
            int_name = splt_ln[1]
            interface_ids = splt_ln[3:]
            interfaces[int_name] = interface_ids
    logging.info(f"Succesfully extracted {len(interfaces)} interface clusters")
    return interfaces


def search_interfaces(interfaces, uniprot_id):
    """
    Searches interfaces for uniprot id.

    Parameters
    ----------
    interfaces : dict
        dictionary of interfaces
    uniprot_id : str
        uniprot id to search for

    Returns
    -------
    found : bool
        True if uniprot id is found

    sth_else : bool
        True if uniprot id is found clustered with something else
    """
    found = False
    sth_else = False

    for int_id in interfaces:
        uniprot_ids = [
            interface.split("-")[0] for interface in interfaces[int_id]
        ]
        unq_uniprot_ids = list(set(uniprot_ids))
        if uniprot_id in unq_uniprot_ids:
            found = True
            logging.info(f"found uniprot id {uniprot_id}")
            if len(unq_uniprot_ids) > 1:
                sth_else = True
                logging.info(
                    f"uniprot id {uniprot_id} is clustered with sth else"
                )
    return found, sth_else


def cycle_run(arctic_io):
    """run arctic for the set of pdb files."""
    stats = {
        # number of times arctic does something
        "n_ids": 0,
        # number of effective runs (no antibodies)
        "n_runs": 0,
        # number of times opponent uniprot id is found
        "n_founds": 0,
        # number of times opponent uniprot id is found clustered
        #  with something else
        "n_clustered_else": 0,
        # number of times arctic3d fails. Ideally zero.
        "n_failures": 0,
    }

    for key in list(arctic_io.keys()):
        rec_pdb = arctic_io[key]["receptor_pdb"]
        self_uni = arctic_io[key]["self_uniprot_id"]
        paired_uni = arctic_io[key]["paired_uniprot_id"]
        logging.info(
            f"processing key {key} (receptor_pdb {rec_pdb}):"
            f" self uni {self_uni} paired uni {paired_uni}"
        )
        # run arctic3d in full mode
        stats["n_ids"] += 1
        pdb_string = ""
        if arctic_io[key]["receptor_pdb"] == "2BBA_A":
            pdb_string += f"--pdb_to_use={rec_pdb[:4]}"
        arctic_folder = f"arctic3d-{self_uni}"
        cmd = (
            "arctic3d"
            f" {self_uni} --full"
            f" --run_dir={arctic_folder} {pdb_string}"
        )
        if Path(arctic_folder).exists():
            logging.info(f"folder {arctic_folder} already exists")
        else:
            output = run_arctic(cmd)
            stats["n_runs"] += 1
            if output == "FAILED":
                stats["n_failures"] += 1
                continue
        int_file = Path(arctic_folder, INT_FILENAME)
        if os.path.exists(int_file):
            interfaces = read_interface_file(int_file)
            arctic_io[key]["arctic_output"] = interfaces
            found, sth_else = search_interfaces(interfaces, paired_uni)
            if found:
                stats["n_founds"] += 1
            if sth_else:
                stats["n_clustered_else"] += 1
            if found:
                # re-running arctic excluding the paired uniprot id
                cmd = (
                    "arctic3d"
                    f" {self_uni} --full"
                    f" --out_partner={paired_uni}"
                    f" --run_dir=excl{paired_uni}-arctic3d-{self_uni}"
                )
                output = run_arctic(cmd)
        else:
            logging.info(
                f"Warning: no interface file found for key {key}:"
                f" receptor pdb {rec_pdb} uniprot {self_uni}"
            )
        logging.info(f"Final stats : {stats}")


# saving ouput
def write_output_file(arctic_io, output_file):
    """Write the output file."""
    logging.info(f"Saving output to {output_file}.")

    with open(output_file, "w") as wfile:
        for key in arctic_io:
            wfile.write(f"{key} {arctic_io[key]} {os.linesep}")


def main():
    # monitoring time
    start_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument(
        "output_folder"
    )  # directory where the benchmark must be built
    args = parser.parse_args()

    bm5_input_file = args.input_file
    output_folder = args.output_folder
    if not os.path.exists(output_folder):
        raise Exception("output folder does not exist")

    output_file = Path(args.output_folder, "bm5_arctic_output.txt")

    arctic_io = get_arctic_io(bm5_input_file)
    os.chdir(output_folder)
    cycle_run(arctic_io)
    write_output_file(arctic_io, output_file)

    elap_time = time.time() - start_time
    logging.info(f"Overall elapsed time: {elap_time:.3f} seconds")


if __name__ == "__main__":
    main()
