"""
Get subcellular localisation of arctic3d data.

Given clustered_interfaces.out input file it iterates over the different partners to detect their
subcellular localisation (as provided by https://www.ebi.ac.uk/proteins/api/proteins)

USAGE::

    arctic3d_resclust ./example/clustered_interfaces.out

"""
import argparse
import logging
import sys
from pathlib import Path
import time
import numpy as np
import matplotlib.pyplot as plt

from arctic3d.functions import make_request
from arctic3d.modules.output import parse_clusters

log = logging.getLogger("arctic3dlog")
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


def main(input_arg):
    """Main function."""
    log.setLevel("INFO")
    start_time = time.time()

    # check input existence
    input_path = Path(input_arg)
    log.info(f"Input file is {input_path}")
    if not input_path.exists():
        raise Exception("non existing input file")


    clustering_dict = parse_clusters(input_path)
    log.info(f"Retrieved clustering_dict with {len(clustering_dict.keys())} clusters.")

    log.info("Retrieving subcellular localisations...")
    locs = {} # partner-specific localisations. Can be empty
    failed_ids = []
    bins = [] # different localisations retrieved
    for cl_id in clustering_dict.keys():
        for partner in clustering_dict[cl_id]:
            uniprot_id = partner.split("-")[0]
            if uniprot_id not in locs.keys() and uniprot_id not in failed_ids:
                uniprot_url = f"{UNIPROT_API_URL}/{uniprot_id}"
                try:
                    prot_data = make_request(uniprot_url, None)
                except Exception as e:
                    log.warning(f"Could not make UNIPROT request for {uniprot_id}, {e}")
                    failed_ids.append(uniprot_id)
                if 'comments' in prot_data.keys():
                    for tup in prot_data['comments']:
                        if tup['type'] == "SUBCELLULAR_LOCATION":
                            if uniprot_id not in locs.keys():
                                locs[uniprot_id] = []
                            for loc in tup['locations']:
                                splt_list = [el.strip() for el in loc['location']['value'].split(",")]
                                locs[uniprot_id].extend(splt_list)
                                for el in splt_list:
                                    if el not in bins:
                                        bins.append(el)
    elap_time = round((time.time() - start_time), 3)
    log.info(f"Subcellular localisation took {elap_time} seconds")
    log.info(f"Retrieved subcellular localisation for {len(locs.keys())} partners.")

    log.info(f"Unique subcellular localisations {bins}")
    
    cl_bins = {}
    # according to the clustering
    for cl_id in clustering_dict.keys():
        cl_bins[cl_id] = {}
        processed_uniprot_ids = []
        for partner in clustering_dict[cl_id]:
            uniprot_id = partner.split("-")[0]
            if uniprot_id in locs.keys() and uniprot_id not in processed_uniprot_ids:
                processed_uniprot_ids.append(uniprot_id)
                for subloc in locs[uniprot_id]:
                    if subloc not in cl_bins[cl_id].keys():
                        cl_bins[cl_id][subloc] = 0
                    cl_bins[cl_id][subloc] += 1
            
    log.info(f"cl_bins {cl_bins}")

    # histograms
    for cluster in cl_bins.keys():
        if cl_bins[cluster] != {}:
            sort_dict = {k: v for k, v in sorted(cl_bins[cluster].items(), reverse=True, key=lambda item: item[1])}
            log.info(f"sort_dict for cluster {cluster} = {sort_dict}")
            plt.figure(figsize=(12,12))
            plt.title(f"cluster_{cluster}")
            plt.bar(sort_dict.keys(), sort_dict.values(), color='g')
            plt.xticks(rotation=45, fontsize=14)
            plt.tight_layout()
            plt.savefig(f"cluster_{cluster}.png")
            plt.close()

    elap_time = round((time.time() - start_time), 3)
    log.info(f"arctic3d-localise run took {elap_time} seconds")


if __name__ == "__main__":
    sys.exit(maincli())
