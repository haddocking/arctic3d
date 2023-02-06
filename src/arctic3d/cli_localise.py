"""
Get subcellular location of arctic3d data.

Given clustered_interfaces.out input file it iterates over the different partners to detect their
subcellular location (as provided by https://www.ebi.ac.uk/proteins/api/proteins)

USAGE::

    arctic3d_localise ./example/clustered_interfaces.out

Use the run_dir parameter if you want to specify a specific output directory::

    arctic3d_localise ./example/clustered_interfaces.out --run_dir=arctic3d-localise-example

Use the out_partner parameter to exclude one or more uniprot IDs from the search::

    arctic3d_localise ./example/clustered_interfaces.out --out_partner=P00760,P00761

It is possible to retrieve information from quickGO instead of using the standard uniprot subcellular location.

QuickGO possesses information location (labelled as C), function (F), and biological process (P)::
    
    arctic3d_localise ./example/clustered_interfaces.out --quickgo=F
"""
import argparse
import logging
import os
import shutil
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt

from arctic3d.functions import make_request
from arctic3d.modules.interface import parse_out_partner
from arctic3d.modules.output import parse_clusters, setup_output_folder, write_dict

LOGNAME = "arctic3d_localise.log"
logging.basicConfig(filename=LOGNAME)
log = logging.getLogger(LOGNAME)
ch = logging.StreamHandler()
formatter = logging.Formatter(
    " [%(asctime)s %(module)s:L%(lineno)d %(levelname)s] %(message)s"
)
ch.setFormatter(formatter)
log.addHandler(ch)

UNIPROT_API_URL = "https://www.ebi.ac.uk/proteins/api/proteins"

argument_parser = argparse.ArgumentParser()
argument_parser.add_argument(
    "input_arg",
    help="Input clustered_interfaces.out file",
)

argument_parser.add_argument(
    "--run_dir",
    help="directory where to store the run",
    default="arctic3d-localise",
)

argument_parser.add_argument(
    "--out_partner",
    help="Set of comma-separated partner IDs to exclude from the search",
)

argument_parser.add_argument(
    "--quickgo",
    help="Use quickgo () information instead of uniprot",
    required=False,
    choices=["C", "F", "P"],
)


def get_uniprot_subcellular_location(prot_data):
    """
    retrieve uniprot subcellular location

    Parameters
    ----------
    prot_data : dict
        uniprot API parsed output dictionary

    Returns
    -------
    locs : list
        list of uniprot subcellular locations
    """
    locs = []
    if "comments" in prot_data.keys():
        for tup in prot_data["comments"]:
            if tup["type"] == "SUBCELLULAR_LOCATION":
                for loc in tup["locations"]:
                    splt_list = [
                        el.strip()
                        for el in loc["location"]["value"].split(",")
                    ]
                    locs.extend(splt_list)
    return locs


def get_quickgo_information(prot_data, quickgo_key):
    """
    retrieve quickgo information

    Parameters
    ----------
    prot_data : dict
        uniprot API parsed output dictionary

    quickgo_key : str
        one among C, F and P

    Returns
    -------
    locs : list
        list of locations (C), functions (F) or biol. processes (B)
    """
    locs = []
    if "dbReferences" in prot_data.keys():
        for tup in prot_data["dbReferences"]:
            if tup["type"] == "GO":
                loc = tup["properties"]["term"]
                if loc.startswith(quickgo_key):
                    locs.append(loc[2:])
    return locs


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


def main(input_arg, run_dir, out_partner, quickgo):
    """Main function."""
    log.setLevel("INFO")
    start_time = time.time()
    if quickgo:
        log.info("Running arctic3d_localise with QUICKGO information")
    else:
        log.info("Running arctic3d_localise with UNIPROT information")

    # check input existence
    input_files = {"cl_filename": Path(input_arg)}
    log.info(f"Input file is {input_files['cl_filename']}")
    if not input_files["cl_filename"].exists():
        raise Exception("non existing input file")

    input_files = setup_output_folder(None, input_files, run_dir)

    # parsing arctic3d clustering output
    clustering_dict = parse_clusters(input_files["cl_filename"])
    log.info(
        "Retrieved clustering_dict with"
        f" {len(clustering_dict.keys())} clusters."
    )

    # parsing out_partner string
    out_partner_set = parse_out_partner(out_partner)
    if out_partner:
        log.info(f"Excluding uniprot IDs {out_partner_set}.")

    log.info("Retrieving subcellular locations...(this may take a while)")
    locs = {}  # partner-specific locations. Can be empty
    failed_ids, none_ids = (
        [],
        [],
    )  # uniprot ids whose call failed/returned None location
    bins = []  # different locations retrieved
    for cl_id in clustering_dict.keys():
        for partner in clustering_dict[cl_id]:
            uniprot_id = partner.split("-")[0]
            # ugly if-clause to avoid calling uniprot with one of the ids to be excluded
            if (
                uniprot_id not in locs.keys()
                and uniprot_id not in failed_ids
                and uniprot_id not in none_ids
                and uniprot_id not in out_partner_set
            ):
                # uniprot call
                log.info(f"calling uniprot with uniprot ID {uniprot_id}")
                uniprot_url = f"{UNIPROT_API_URL}/{uniprot_id}"
                try:
                    prot_data = make_request(uniprot_url, None)
                except Exception as e:
                    log.warning(
                        f"Could not make UNIPROT request for {uniprot_id}, {e}"
                    )
                    failed_ids.append(uniprot_id)
                    continue
                # parsing
                if quickgo:
                    locations = get_quickgo_information(
                        prot_data, quickgo_key=quickgo
                    )
                else:
                    locations = get_uniprot_subcellular_location(prot_data)

                if locations == []:
                    log.info(f"no location retrieved for {uniprot_id}")
                    none_ids.append(uniprot_id)
                else:
                    log.info(f"location retrieved for {uniprot_id}")
                    locs[uniprot_id] = locations
                    # append to bins
                    for location in locations:
                        if location not in bins:
                            bins.append(location)
    elap_time = round((time.time() - start_time), 3)
    log.info(f"Subcellular location retrieval took {elap_time} seconds")
    log.info(f"{len(failed_ids)} partners failed uniprot calls.")
    log.info(f"{len(none_ids)} contain None subcellular location information.")
    log.info(
        f"Retrieved subcellular location for {len(locs.keys())} partners."
    )

    log.info(f"Unique subcellular locations {bins}")

    # writing locations to file
    loc_filename = "Subcellular_locations.txt"
    log.info(f"Saving retrieved subcellular locations to file {loc_filename}")
    write_dict(locs, loc_filename, keyword="Subcellular location", sep=",")

    # creating the histograms according to the clustering
    cl_bins = {}
    for cl_id in clustering_dict.keys():
        cl_bins[cl_id] = {}
        processed_uniprot_ids = []
        for partner in clustering_dict[cl_id]:
            uniprot_id = partner.split("-")[0]
            if (
                uniprot_id in locs.keys()
                and uniprot_id not in processed_uniprot_ids
            ):
                processed_uniprot_ids.append(uniprot_id)
                for subloc in locs[uniprot_id]:
                    if subloc not in cl_bins[cl_id].keys():
                        cl_bins[cl_id][subloc] = 0
                    cl_bins[cl_id][subloc] += 1

    log.info("Plotting cluster locations...")
    # plotting histograms
    for cluster in cl_bins.keys():
        if cl_bins[cluster] != {}:
            sort_dict = {
                k: v
                for k, v in sorted(
                    cl_bins[cluster].items(), key=lambda item: item[1]
                )
            }
            labels = list(sort_dict.keys())
            values = list(sort_dict.values())
            gap = (max(values) - min(values)) // 12 + 1
            xints = range(min(values), max(values) + 1, gap)
            plt.figure(figsize=(12, 12))
            plt.title(f"cluster {cluster}", fontsize=24)
            plt.barh(labels, values, height=0.3, color="g")
            plt.xticks(xints, fontsize=18)
            plt.yticks(fontsize=14)
            plt.tight_layout()
            fig_fname = f"cluster_{cluster}.png"
            plt.savefig(fig_fname)
            log.info(f"Figure {fig_fname} created")
            plt.close()

    # saving histograms
    os.mkdir("histograms")
    for cl_id in cl_bins.keys():
        if cl_bins[cl_id] != {}:
            log.info(f"writing histogram for cluster {cl_id}")
            sort_dict = {
                k: v
                for k, v in sorted(
                    cl_bins[cl_id].items(),
                    key=lambda item: item[1],
                    reverse=True,
                )
            }
            histo_file = Path("histograms", f"cluster_{cl_id}.tsv")
            with open(histo_file, "w") as wfile:
                labels = list(sort_dict.keys())
                values = list(sort_dict.values())
                for n in range(len(labels)):
                    wfile.write(f"{labels[n]}\t{values[n]}{os.linesep}")
        else:
            log.warning(f"cluster {cl_id} empty: will be discarded")

    elap_time = round((time.time() - start_time), 3)
    log.info(f"arctic3d_localise run took {elap_time} seconds")
    shutil.move(f"../{LOGNAME}", LOGNAME)


if __name__ == "__main__":
    sys.exit(maincli())
