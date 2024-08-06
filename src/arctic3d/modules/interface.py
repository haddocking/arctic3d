"Interface library."
import logging
import os

import jsonpickle

from arctic3d.functions import make_request
from arctic3d.modules.url import INTERFACE_URL, LIGAND_URL

log = logging.getLogger("arctic3d.log")

# maximum number of interfaces in interface file
MAX_INTERFACES = 10000


def parse_out_partner(out_partner_string):
    """
    Parse the string of uniprot IDs to exclude.

    Parameters
    ----------
    out_partner_string : str
        comma-separated uniprot IDs to exclude.

    Returns
    -------
    out_partner_list : set
        set of uniprot IDs to exclude.
    """
    out_partner_list = []
    if out_partner_string:
        uniprot_split = out_partner_string.split(",")
        for uni in uniprot_split:
            if not uni.isalnum():
                raise Exception(f"Invalid Uniprot ID {uni} in --out_pdb.")
            if uni in out_partner_list:
                log.warning(
                    f"Warning: duplicated pdb entry {uni} in --out_partner."
                )
            else:
                out_partner_list.append(uni)
    return set(out_partner_list)


def parse_out_pdb(out_pdb_string):
    """
    Parse the string of PDB IDs to exclude.

    Parameters
    ----------
    out_pdb_string : str
        comma-separated PDB IDs to exclude.

    Returns
    -------
    out_pdb_set : set
        set of PDB IDs to exclude.
    """
    out_pdb_list = []
    if out_pdb_string:
        pdb_split = out_pdb_string.split(",")
        for pdb in pdb_split:
            if len(pdb) != 4 or not pdb.isalnum():
                raise Exception(f"Invalid PDB ID {pdb} in --out_pdb.")
            lower_pdb = pdb.lower()
            if lower_pdb in out_pdb_list:
                log.warning(
                    f"Warning: duplicated pdb entry {lower_pdb} in --out_pdb."
                )
            else:
                out_pdb_list.append(lower_pdb)
    return set(out_pdb_list)


def parse_interface_line(int_line, ln_num):
    """
    Parses the input interface line according to the following format:

    int_name 1,2,3,6,7

    Parameters
    ----------
    int_line : str
        interface_line
    ln_num : int
        line number

    Returns
    -------
    int_name : str
        name of the interface
    residue_list : list
        list of residues
    """
    splt_ln = int_line.strip().split()
    int_name = splt_ln[0]
    # checking malformed interface
    if len(splt_ln) != 2:
        raise Exception(
            f"Found uncompatible interface at line {ln_num} in interface_file."
        )
    residues_str_list = splt_ln[1].split(",")
    residues_int_list = []
    # checking they are all integers
    for resid_string in residues_str_list:
        if resid_string.isdigit():
            residues_int_list.append(int(resid_string))
        else:
            raise Exception(
                f"Malformed residue {resid_string} at line {ln_num} in"
                " interface_file."
            )
    return int_name, residues_int_list


def read_interface_residues(interface_file):
    """
    Parameters
    ----------
    interface_file : str or Path
        path to the interface file

    Example_file :
    int_1 1,2,3
    int_2 1,2,4,5,6
    ...
    interface_name_N 35,78,18

    Returns
    -------
    interface_residues : dict
        Interface residues.
    """
    interface_dict = {}
    if os.path.exists(interface_file):
        with open(interface_file, "r") as ifile:
            ln_num = 0  # keep track of line number
            for ln in ifile:
                ln_num += 1
                if ln != os.linesep:
                    int_name, residue_list = parse_interface_line(ln, ln_num)
                    interface_dict[int_name] = residue_list
        if ln_num > MAX_INTERFACES:
            raise Exception(
                f"Number of interfaces ({ln_num}) higher than threshold"
                f" ({MAX_INTERFACES})."
            )
    else:
        raise Exception(f"interface_file {interface_file} does not exist")
    return interface_dict


def get_interface_residues(
    uniprot_id,
    out_partner_string,
    out_pdb_string,
    full,
    ligand,
    interface_data=None,
):
    """
    Get interface residues.

    Parameters
    ----------
    uniprot_id : str
        Uniprot ID.
    out_partner_string : str or None
        comma-separated uniprot IDs to exclude.
    out_pdb_string : str or None
        comma-separated pdb IDs to exclude.
    full : bool
        consider full information in interface retrieval
    interface_data : str or Path or None
        interface data .json file
    ligand : str
        retrieve ligand binding residues

    Returns
    -------
    interface_residues : dict
        Interface residues.

    """
    # get the uniprot IDs to exclude
    out_partner_set = parse_out_partner(out_partner_string)
    # get the PDB IDs to exclude
    out_pdb_set = parse_out_pdb(out_pdb_string)

    interface_dict = {}
    if interface_data:
        if ligand == "no":
            try:
                interface_api_data = jsonpickle.decode(
                    open(interface_data, "r").read()
                )
            except Exception as e:
                log.warning(
                    "Could not read input interface_data"
                    f" {interface_data}, {e}"
                )
                return interface_dict
        elif ligand == "yes":
            try:
                interface_lig_api_data = jsonpickle.decode(
                    open(interface_data, "r").read()
                )
            except Exception as e:
                log.warning(
                    "Could not read input interface_data"
                    f" {interface_data}, {e}"
                )
                return interface_dict
    else:
        if ligand in ["no", "both"]:
            url = f"{INTERFACE_URL}/{uniprot_id}"
            try:
                interface_api_data = make_request(url, None)
            except Exception as e:
                log.warning(
                    "Could not make InterfaceResidues request for"
                    f" {uniprot_id}, {e}"
                )
                return interface_dict
        if ligand in ["yes", "both"]:
            url = f"{LIGAND_URL}/{uniprot_id}"
            try:
                interface_lig_api_data = make_request(url, None)
            except Exception as e:
                log.warning(
                    f"Could not make LigandSites request for {uniprot_id}, {e}"
                )
                return interface_dict

    if ligand in ["no", "both"]:
        if interface_api_data and len(interface_api_data) != 0:
            interface_dict = parse_interface_data(
                uniprot_id=uniprot_id,
                interface_data=interface_api_data,
                out_partner_set=out_partner_set,
                out_pdb_set=out_pdb_set,
                full=full,
            )
    if ligand in ["yes", "both"]:
        if interface_lig_api_data and len(interface_lig_api_data) != 0:
            interface_lig_dict = parse_interface_data(
                uniprot_id=uniprot_id,
                interface_data=interface_lig_api_data,
                out_partner_set=out_partner_set,
                out_pdb_set=out_pdb_set,
                full=full,
            )
    if ligand == "yes":
        interface_dict = interface_lig_dict
    if ligand == "both":
        interface_dict.update(interface_lig_dict)

    return interface_dict


def parse_interface_data(
    uniprot_id, interface_data, out_partner_set, out_pdb_set, full
):
    """
    Parse interface data.

    Parameters
    ----------
    uniprot_id : str
        Uniprot ID.
    interface_data : dict
        Interface data.
    out_partner_list : set
        Set of Uniprot IDs to exclude. Empty if none.
    out_pdb_set : set
        Set of PDB files to exclude.
    full : bool
        Consider full information in interface retrieval.

    Returns
    -------
    interface_dict : dict
        Interface residue dictionary.
        *example* : {
            partner_uniprotid_1: [1,2,3],
            partner_uniprotid_2: [20,22,23]
            }
    """
    interface_dict = {}
    for element in interface_data[uniprot_id]["data"]:
        partner_uniprotid = element["accession"]
        if partner_uniprotid not in out_partner_set:
            # log.info(f"Parsing partner uniprot ID {partner_uniprotid}")
            for residue_entry in element["residues"]:
                start = residue_entry["startIndex"]
                end = residue_entry["endIndex"]
                if full is False:
                    # consider all pdb files as a whole
                    accept = True
                    key = partner_uniprotid
                    if out_pdb_set:
                        int_pdbs_list = residue_entry["interactingPDBEntries"]
                        int_pdbs = set(
                            [data["pdbId"] for data in int_pdbs_list]
                        )
                        if int_pdbs.issubset(out_pdb_set):
                            # all interacting pdbs must be excluded
                            accept = False
                    if accept:
                        if key not in interface_dict.keys():
                            interface_dict[key] = []
                        for interface_res in range(start, end + 1):
                            interface_dict[partner_uniprotid].append(
                                interface_res
                            )
                else:
                    # iterate over pdb records
                    for pdb_record in residue_entry["interactingPDBEntries"]:
                        # entries can have missing chainIds field, especially
                        #   for PRD_* like uniprot IDs
                        chain_ids = [""]
                        if "chainIds" in pdb_record.keys():
                            chain_ids = pdb_record["chainIds"].split(",")
                        # if there are two or more chainIds, we discard the
                        #   current entry
                        if pdb_record["pdbId"] not in out_pdb_set:
                            for chain_id in chain_ids:
                                key = (
                                    f"{partner_uniprotid}"
                                    f"-{pdb_record['pdbId']}-{chain_id}"
                                )
                                if key not in interface_dict.keys():
                                    interface_dict[key] = []
                                for interface_res in range(start, end + 1):
                                    interface_dict[key].append(interface_res)
        else:
            log.info(
                f"found uniprot ID {partner_uniprotid}. It will be discarded."
            )

    log.info(f"{len(interface_dict.keys())} retrieved interfaces.")
    return interface_dict
