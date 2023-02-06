import argparse
import logging
import os
import shlex
import shutil
import subprocess
import time
from pathlib import Path

import pandas as pd

INT_FILENAME = "clustered_interfaces.out"
RES_FILENAME = "clustered_residues.out"


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

        arctic_io[receptor_pdb] = {
            "complex_pdb": complex_pdb,
            "self_uniprot_id": receptor_uniprot,
            "paired_uniprot_id": ligand_uniprot,
            "arctic_output": {},
        }

        arctic_io[ligand_pdb] = {
            "complex_pdb": complex_pdb,
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
    # out = p.stdout.decode("utf-8").split(os.linesep)
    # check if failed
    if p.returncode != 0:
        print("warning: ARCTIC failed")
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
            print(f"interface ids {interface_ids}")
            interfaces[int_name] = interface_ids
    return interfaces


def search_interfaces(interfaces, uniprot_id):
    found = False
    sth_else = False
    for int_id in interfaces:
        if uniprot_id in interfaces[int_id]:
            found = True
            print(f"found uniprot id {uniprot_id}")
            if len(interfaces[int_id]) > 1:
                sth_else = True
                print(f"uniprot id {uniprot_id} is clustered with sth else")
    return found, sth_else


def cycle_run(arctic_io):
    """run arctic for the set of pdb files."""
    stats = {
        # number of times arctic runs (no antibodies)
        "n_runs": 0,
        # number of times opponent uniprot id is found
        "n_founds": 0,
        # number of times opponent uniprot id is found clustered
        #  with something else
        "n_clustered_else": 0,
        # number of times arctic3d fails. Ideally zero.
        "n_failures": 0,
    }

    for pdb in list(arctic_io.keys()):
        print(f"processing pdb {pdb}")
        complex_code = arctic_io[pdb]["complex_pdb"][:4]
        # run arctic3d in full mode
        stats["n_runs"] += 1
        cmd = f"arctic3d {arctic_io[pdb]['self_uniprot_id']}"
        output = run_arctic(cmd)
        if output == "FAILED":
            stats["n_failures"] += 1
        else:
            if os.path.exists(INT_FILENAME):
                interfaces = read_interface_file(INT_FILENAME)
                arctic_io[pdb]["arctic_output"] = interfaces
                found, sth_else = search_interfaces(
                    interfaces, arctic_io[pdb]["paired_uniprot_id"]
                )
                if found:
                    stats["n_founds"] += 1
                if sth_else:
                    stats["n_clustered_else"] += 1
                # mkdir if not already present
                if not os.path.exists(complex_code):
                    os.mkdir(complex_code)
                # copy stuff to complex_code directory
                int_file = Path(
                    complex_code,
                    f"clustered_interfaces_{pdb[:4]}"
                    f"_{arctic_io[pdb]['self_uniprot_id']}_full.out",
                )
                res_file = Path(
                    complex_code,
                    f"clustered_residues_{pdb[:4]}"
                    f"_{arctic_io[pdb]['self_uniprot_id']}_full.out",
                )
                shutil.copy(INT_FILENAME, int_file)
                shutil.copy(RES_FILENAME, res_file)
                os.unlink(RES_FILENAME)
                os.unlink(INT_FILENAME)
                if found:
                    # re-running arctic excluding the paired uniprot id
                    cmd = (
                        "arctic3d"
                        f" {arctic_io[pdb]['self_uniprot_id']} --out_uniprot"
                        f" {arctic_io[pdb]['paired_uniprot_id']}"
                    )
                    output = run_arctic(cmd)
                    if output == "SUCCESS" and os.path.exists(INT_FILENAME):
                        int_file = Path(
                            complex_code,
                            f"clustered_interfaces_{pdb[:4]}"
                            f"_{arctic_io[pdb]['self_uniprot_id']}_excl.out",
                        )
                        res_file = Path(
                            complex_code,
                            f"clustered_residues_{pdb[:4]}"
                            f"_{arctic_io[pdb]['self_uniprot_id']}_excl.out",
                        )
                        shutil.copy(INT_FILENAME, int_file)
                        shutil.copy(RES_FILENAME, res_file)
                        os.unlink(RES_FILENAME)
                        os.unlink(INT_FILENAME)
            else:
                print(
                    f"Warning: no interface file found for pdb {pdb} uniprot"
                    f" {arctic_io[pdb]['self_uniprot_id']}"
                )
        print(f"Final stats : {stats}")


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
    print(f"Overall elapsed time: {elap_time:.3f} seconds")


if __name__ == "__main__":
    main()
