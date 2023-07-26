"""
Get subcellular location of arctic3d data.

Given clustered_interfaces.out input file it iterates over the
    different partners to detect their
subcellular location (as provided by
    https://www.ebi.ac.uk/proteins/api/proteins)

USAGE::

    arctic3d-localise ./example/clustered_interfaces.out

Use the run_dir parameter if you want to specify a specific output directory::

    arctic3d-localise ./example/clustered_interfaces.out \
        --run_dir=arctic3d-localise-example

Use the out_partner parameter to exclude one or more uniprot IDs
from the search::

    arctic3d-localise ./example/clustered_interfaces.out \
        --out_partner=P00760,P00761

It is possible to retrieve information from quickGO instead of using
the standard uniprot subcellular location.

QuickGO possesses information location (labelled as C), function (F),
and biological process (P)::

    arctic3d-localise ./example/clustered_interfaces.out \
        --quickgo=F

By default the CLI will output html files and png images. It is possible
to specify the output format with the format parameter::

    arctic3d-localise ./example/clustered_interfaces.out \
        --format=pdf

The scale parameter can be used to scale the images::

    arctic3d-localise ./example/clustered_interfaces.out \
        --scale=2.0
"""
import argparse
import os
import sys
import time
from pathlib import Path

from arctic3d import log

from arctic3d.functions import make_request
from arctic3d.modules.interface import parse_out_partner
from arctic3d.modules.log import add_log_for_CLI
from arctic3d.modules.output import (
    #    create_barplot,
    create_barplotly,
    create_output_folder,
    parse_clusters,
    setup_output_folder,
    write_dict,
)

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
    help="Use quickgo (https://www.ebi.ac.uk/QuickGO/annotations) information "
    "instead of uniprot",
    required=False,
    choices=["C", "F", "P"],
)

argument_parser.add_argument(
    "--weight",
    help="Weight histograms according to uniprot. The more populated an"
    " UNIPROT record, the less weight it will have in the histogram.",
    required=False,
    choices=["yes", "no"],
    default="no",
)

argument_parser.add_argument(
    "--format",
    help="produce images in the desired format",
    required=False,
    type=str,
    default="png",
    choices=["png", "pdf", "svg", "jpeg", "webp"],
)

argument_parser.add_argument(
    "--scale", help="scale for images", required=False, type=float, default=4.0
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


def get_sorted_dict(cluster_bins, reverse=False):
    """
    Sorts dictionary

    Parameters
    ----------
    cluster_bins : dict
        dictionary to sort according to its values

    reverse : bool
        use reverse order

    Returns
    -------
    sort_dict : dict
        sorted dictionary
    """
    sort_dict = {
        k: v
        for k, v in sorted(
            cluster_bins.items(), key=lambda item: item[1], reverse=reverse
        )
    }
    return sort_dict


def create_histograms(uniprot_clustering_dict, locs, weight):
    """
    Create histograms.

    Parameters
    ----------
    uniprot_clustering_dict : dict
        uniprot-based clustering dictionary
    locs : dict
        locations dictionary, with uniprot IDs as keys
    weight : str
        weight histograms according to uniprot.
    """
    cl_bins = {}
    for cl_id in uniprot_clustering_dict.keys():
        cl_bins[cl_id] = {}
        for uniprot_id in uniprot_clustering_dict[cl_id]:
            if uniprot_id in locs.keys():
                # adjusting weight
                if weight == "yes":
                    bin_weight = 1 / len(locs[uniprot_id])
                else:
                    bin_weight = 1
                # looping over locations
                for subloc in locs[uniprot_id]:
                    if subloc not in cl_bins[cl_id].keys():
                        cl_bins[cl_id][subloc] = 0
                    cl_bins[cl_id][subloc] += bin_weight
    return cl_bins


def get_uniprot_dict(clustering_dict, out_partner_set):
    """
    Get uniprot dictionary
    """
    uniprot_clustering_dict = {}
    uniprot_set = []
    for cl_id in clustering_dict.keys():
        uniprot_clustering_dict[cl_id] = []
        # loop over partners
        for partner in clustering_dict[cl_id]:
            uniprot_id = partner.split("-")[0]
            if uniprot_id not in out_partner_set:
                if uniprot_id not in uniprot_clustering_dict[cl_id]:
                    uniprot_clustering_dict[cl_id].append(uniprot_id)
                uniprot_set.append(uniprot_id)
    uniprot_set = list(set(uniprot_set))
    return uniprot_clustering_dict, uniprot_set


def call_uniprot(uniprot_id):
    """
    Call uniprot API

    Parameters
    ----------
    uniprot_id : str
        uniprot ID

    Returns
    -------
    prot_data : dict
        uniprot API parsed output dictionary
    """
    # uniprot call
    log.info(f"calling uniprot with uniprot ID {uniprot_id}")
    uniprot_url = f"{UNIPROT_API_URL}/{uniprot_id}"
    try:
        prot_data = make_request(uniprot_url, None)
    except Exception as e:
        log.warning(f"Could not make UNIPROT request for {uniprot_id}, {e}")
        prot_data = None
    return prot_data


def get_locations(prot_data, quickgo):
    """
    Get locations from uniprot or quickgo

    Parameters
    ----------
    prot_data : dict
        uniprot API parsed output dictionary
    quickgo : str
        one among C, F and P

    Returns
    -------
    locations : list
    """
    if quickgo:
        locations = get_quickgo_information(prot_data, quickgo_key=quickgo)
    else:
        locations = get_uniprot_subcellular_location(prot_data)
    return locations


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
    run_dir,
    out_partner,
    quickgo,
    weight,
    format,
    scale,
    log_level="DEBUG",
):
    """Main function."""
    log.setLevel(log_level)
    start_time = time.time()

    # create output folder
    run_dir_path = create_output_folder(run_dir)
    # logging
    log_file = Path(run_dir_path, "arctic3d-localise.log")
    add_log_for_CLI(log, log_level, log_file)

    # property name
    prop_name = "location"
    if quickgo:
        if quickgo == "F":
            prop_name = "function"
        elif quickgo == "P":
            prop_name = "biological process"
        log.info(
            f"Running arctic3d-localise with QUICKGO {prop_name} information"
        )
    else:
        log.info(
            f"Running arctic3d-localise with UNIPROT {prop_name} information"
        )

    # check input existence
    input_files = {"cl_filename": Path(input_arg)}
    log.info(f"Input file is {input_files['cl_filename']}")
    if not input_files["cl_filename"].exists():
        raise Exception("non existing input file")

    input_files = setup_output_folder(run_dir_path, input_files)

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

    # writing down uniprot-based clustering (to make sense of the data)
    uniprot_clustering_dict, uniprot_set = get_uniprot_dict(
        clustering_dict, out_partner_set
    )

    # write down uniprot_clustering_dict
    write_dict(
        uniprot_clustering_dict,
        "uniprot_clustering.out",
        keyword="Cluster",
        sep=",",
    )

    log.info(
        f"Retrieving {prop_name} of {len(uniprot_set)} partners..."
        "(this may take a while)"
    )
    locs = {}  # partner-specific locations. Can be empty
    failed_ids, none_ids = (
        [],
        [],
    )  # uniprot ids whose call failed/returned None location
    bins = []  # different locations retrieved
    for uniprot_id in uniprot_set:
        prot_data = call_uniprot(uniprot_id)
        if prot_data is None:
            failed_ids.append(uniprot_id)
            continue
        # parsing
        locations = get_locations(prot_data, quickgo)

        if locations == []:
            log.info(f"no {prop_name} retrieved for {uniprot_id}")
            none_ids.append(uniprot_id)
        else:
            log.info(f"{prop_name} retrieved for {uniprot_id}")
            locs[uniprot_id] = locations
            # append to bins
            for location in locations:
                if location not in bins:
                    bins.append(location)

    elap_time = round((time.time() - start_time), 3)
    log.info(f"{prop_name} retrieval took {elap_time} seconds")
    log.info(f"{len(failed_ids)} partners failed uniprot calls.")
    log.info(f"{len(none_ids)} contain None {prop_name} information.")
    log.info(f"Retrieved {prop_name} for {len(locs.keys())} partners.")

    log.info(f"Unique {prop_name} entries: {bins}")

    # writing locations to file
    loc_filename = f"{prop_name}.txt"
    log.info(f"Saving retrieved {prop_name} to file {loc_filename}")
    write_dict(locs, loc_filename, keyword=prop_name, sep=",")

    # creating the histograms according to the clustering
    cl_bins = create_histograms(uniprot_clustering_dict, locs, weight)

    log.info("Plotting cluster locations...")
    # plotting histograms
    for cluster in cl_bins.keys():
        if cl_bins[cluster] != {}:
            # get sorted dictionary
            sort_dict = get_sorted_dict(cl_bins[cluster])
            # plot
            # create_barplot(cluster, sort_dict, max_labels=50)
            create_barplotly(cluster, sort_dict, format, scale, max_labels=25)

    # saving histograms
    os.mkdir("histograms")
    for cl_id in cl_bins.keys():
        if cl_bins[cl_id] != {}:
            log.info(f"writing histogram for cluster {cl_id}")
            sort_dict = get_sorted_dict(cl_bins[cl_id], reverse=True)
            histo_file = Path("histograms", f"cluster_{cl_id}.tsv")
            with open(histo_file, "w") as wfile:
                labels = list(sort_dict.keys())
                values = list(sort_dict.values())
                for n in range(len(labels)):
                    wfile.write(f"{labels[n]}\t{values[n]}{os.linesep}")
        else:
            log.warning(f"cluster {cl_id} empty: will be discarded")

    elap_time = round((time.time() - start_time), 3)
    log.info(f"arctic3d-localise run took {elap_time} seconds")


if __name__ == "__main__":
    sys.exit(maincli())
