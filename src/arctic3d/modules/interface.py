import logging
import os

import jsonpickle

from arctic3d.functions import make_request

log = logging.getLogger("arctic3dlog")

INTERFACE_URL = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot/interface_residues"
LIGAND_URL = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot/ligand_sites/"
# maximum number of interfaces in interface file
MAX_INTERFACES = 10000
# ALLPDB_URL = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot"
# QUERY_DB = {}  # this probably can be removed
# PROTEINS_URL = "https://www.ebi.ac.uk/proteins/api/proteins"
# Ig_words = [
#     "Ig",
#     "IG",
#     "Immunoglobulin",
#     "immunoglobulin",
#     "IMMUNOGLOBULIN",
#     "Antibody",
#     "antibody",
#     "ANTIBODY",
#     "SAB",
# ]
# filters_list = ["", "", "", "", "", ""]
# NEGATOME_URL = (
#     "http://mips.helmholtz-muenchen.de/proj/ppi/negatome/combined_stringent.txt"
# )
# PDB_MOLECULES_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules"
# COV_CUTOFF = 0.8
# PDB_DB = []


def parse_out_uniprot(out_uniprot_string):
    """
    Parse the string of uniprot IDs to exclude.

    Parameters
    ----------
    out_uniprot_string : str
        comma-separated uniprot IDs to exclude.

    Returns
    -------
    out_uniprot_list : set
        set of uniprot IDs to exclude.
    """
    out_uniprot_list = []
    if out_uniprot_string:
        uniprot_split = out_uniprot_string.split(",")
        for uni in uniprot_split:
            if not uni.isalnum():
                raise Exception(f"Invalid Uniprot ID {uni} in --out_pdb.")
            if uni in out_uniprot_list:
                log.warning(f"Warning: duplicated pdb entry {uni} in --out_uniprot.")
            else:
                out_uniprot_list.append(uni)
    return set(out_uniprot_list)


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
                log.warning(f"Warning: duplicated pdb entry {lower_pdb} in --out_pdb.")
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
                f"Malformed residue {resid_string} at line {ln_num} in interface_file."
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
                f"Number of interfaces ({ln_num}) higher than threshold ({MAX_INTERFACES})."
            )
    else:
        raise Exception(f"interface_file {interface_file} does not exist")
    return interface_dict


def get_interface_residues(
    uniprot_id, out_uniprot_string, out_pdb_string, full, ligand, interface_data=None
):
    """
    Get interface residues.

    Parameters
    ----------
    uniprot_id : str
        Uniprot ID.
    out_uniprot_string : str or None
        comma-separated uniprot IDs to exclude.
    out_pdb_string : str or None
        comma-separated pdb IDs to exclude.
    full : bool
        consider full information in interface retrieval
    interface_data : str or Path or None
        interface data .json file
    ligand : True, False or "both"
        retrieve ligand binding residues 

    Returns
    -------
    interface_residues : dict
        Interface residues.

    """
    # get the uniprot IDs to exclude
    out_uniprot_set = parse_out_uniprot(out_uniprot_string)
    # get the PDB IDs to exclude
    out_pdb_set = parse_out_pdb(out_pdb_string)

    interface_dict = {}
    if not interface_data:
        if ligand in ["no", "both"]:
            url = f"{INTERFACE_URL}/{uniprot_id}"
            try:
                interface_api_data = make_request(url, None)
            except Exception as e:
                log.warning(
                    f"Could not make InterfaceResidues request for {uniprot_id}, {e}"
                )
                return interface_dict
        if ligand in ["yes", "both"]:
            url = f"{LIGAND_URL}/{uniprot_id}"
            try:
                interface_lig_api_data = make_request(url, None)
            except Exception as e:
                log.warning(
                    f"Could not make InterfaceResidues request for {uniprot_id}, {e}"
                )
                return interface_dict
    else:
        try:
            interface_api_data = jsonpickle.decode(open(interface_data, "r").read())
        except Exception as e:
            log.warning(f"Could not read input interface_data {interface_data}, {e}")
            return interface_dict

    if ligand in ["no", "both"]:
        if interface_api_data and len(interface_api_data) != 0:
            interface_dict = parse_interface_data(
                uniprot_id, interface_api_data, out_uniprot_set, out_pdb_set, full
            )
    if ligand in ["yes", "both"]:
        if interface_lig_api_data and len(interface_lig_api_data) != 0:
            interface_lig_dict = parse_interface_data(
                uniprot_id, interface_lig_api_data, [], [], full=full
                #uniprot_id, interface_lig_api_data, out_uniprot_set, out_pdb_set, full=full
            )
    if ligand == "yes":
        interface_dict = interface_lig_dict
    if ligand == "both":
        interface_dict.update(interface_lig_dict)
    
    log.info(f"interface_dict {interface_dict}")

    return interface_dict


def parse_interface_data(
    uniprot_id, interface_data, out_uniprot_set, out_pdb_set, full
):
    """
    Parse interface data.

    Parameters
    ----------
    uniprot_id : str
        Uniprot ID.
    interface_data : dict
        Interface data.
    out_uniprot_list : set
        Set of Uniprot IDs to exclude. Empty if none.
    out_pdb_set : set
        Set of PDB files to exclude.
    full : bool
        Consider full information in interface retrieval.

    Returns
    -------
    interface_dict : dict
        Interface residue dictionary.
            {partner_uniprotid_1: [1,2,3], partner_uniprotid_2: [20,22,23]}
    """
    interface_dict = {}
    for element in interface_data[uniprot_id]["data"]:
        partner_uniprotid = element["accession"]
        if partner_uniprotid not in out_uniprot_set:
            log.info(f"Parsing partner uniprot ID {partner_uniprotid}")
            for residue_entry in element["residues"]:
                start = residue_entry["startIndex"]
                end = residue_entry["endIndex"]
                if full is False:
                    # consider all pdb files as a whole
                    accept = True
                    key = partner_uniprotid
                    if out_pdb_set:
                        int_pdbs_list = residue_entry["interactingPDBEntries"]
                        int_pdbs = set([data["pdbId"] for data in int_pdbs_list])
                        if int_pdbs.issubset(out_pdb_set):
                            # all interacting pdbs must be excluded
                            accept = False
                    if accept:
                        if key not in interface_dict.keys():
                            interface_dict[key] = []
                        for interface_res in range(start, end + 1):
                            interface_dict[partner_uniprotid].append(interface_res)
                else:
                    # iterate over pdb records
                    for pdb_record in residue_entry["interactingPDBEntries"]:
                        # entries can have missing chainIds field, especially for PRD_* like uniprot IDs
                        chain_ids = [""]
                        if "chainIds" in pdb_record.keys():
                            chain_ids = pdb_record["chainIds"].split(",")
                        # if there are two or more chainIds, we discard the current entry
                        if pdb_record["pdbId"] not in out_pdb_set:
                            for chain_id in chain_ids:
                                key = f"{partner_uniprotid}-{pdb_record['pdbId']}-{chain_id}"
                                if key not in interface_dict.keys():
                                    interface_dict[key] = []
                                for interface_res in range(start, end + 1):
                                    interface_dict[key].append(interface_res)
        else:
            log.info(f"found uniprot ID {partner_uniprotid}. It will be discarded.")

    log.info(f"{len(interface_dict.keys())} retrieved interfaces.")
    return interface_dict


# def get_uniprot_length(uniprot_id):
#     """
#     gets the length of the uniprot id sequence

#     Parameters
#     ----------
#     uniprot_id : str
#         uniprot ID

#     Returns
#     -------
#     id_length : int
#         length of uniprot ID
#     """
#     lig_len = None
#     lig_interface_url = f"{INTERFACE_URL}/{uniprot_id}"
#     try:
#         lig_interface_data = make_request(lig_interface_url, None)
#     except Exception as e:
#         log.warning(f"Could not make INTERFACE request for {uniprot_id}, {e}")

#     if lig_interface_data:
#         if len(lig_interface_data) != 0:  # Check if the entry is empty
#             if len(lig_interface_data[uniprot_id]) == 0:  # Check if the entry is empty.
#                 log.warning(f"InterfaceResidues for {uniprot_id} was empty")
#             else:  # There is information
#                 lig_len = len(
#                     lig_interface_data[uniprot_id]["sequence"]
#                 )  # Lenght of the receptor according to Uniprot
#     return lig_len


# def calculate_coverage(pdb_structure, pdb_id, reference_length):
#     """check coverage"""
#     coverage = False
#     log.info(f"checking coverage of pdb id {pdb_id}")
#     if len(pdb_structure["segments"]) == 1:  # *Jesus: This can be improved I guess
#         pdb_len = len(pdb_structure["segments"][0]["pdb_sequence"])
#         pdb_coverage = float(pdb_len / reference_length)
#         log.info(f"coverage is {pdb_coverage}")
#         if pdb_coverage >= COV_CUTOFF:  # Successful receptor coverage
#             coverage = True
#     return coverage


# def check_partners_coverage(receptor_id, rec_len, ligand_id, int_pdbs):
#     """
#     Check if the coverage of the interacting partners surpass the established threshold.

#     Parameters
#     ----------

#     receptor_id : str
#         uniprot ID

#     rec_len : int
#         receptor length

#     ligand_id : str
#         uniprot id of the ligand

#     int_pdbs : list
#         list of interacting pdb ids

#     Returns
#     -------
#     confirmed_pdbs : list
#         list of pdbs that passed the sequence coverage assessment

#     """
#     log.info(f"int_pdbs {int_pdbs}")
#     confirmed_pdbs = None
#     # Get the Uniprot length for the ligand
#     lig_len = get_uniprot_length(ligand_id)
#     if lig_len is None:
#         log.warning("Cannot compute ligand length")
#     else:
#         confirmed_pdbs = []

#         # Information from the PDBs
#         rec_url = f"{ALLPDB_URL}/{receptor_id}"
#         lig_url = f"{ALLPDB_URL}/{ligand_id}"

#         # Check receptor coverage
#         try:
#             rec_allpdb_data = make_request(rec_url, None)
#         except Exception as e:
#             log.warning(f"Could not make AllPDB request for {receptor_id}, {e}")

#         if rec_allpdb_data:
#             rec_pdbs = []
#             # Iterate every pdb for the receptor to check if it is one of the preselected ones
#             for rec_pdb_structure in rec_allpdb_data[receptor_id]["mappings"]:
#                 pdb_id = rec_pdb_structure["entry_id"]
#                 if pdb_id in int_pdbs:
#                     cov = calculate_coverage(rec_pdb_structure, pdb_id, rec_len)
#                     if cov >= COV_CUTOFF:  # Successful receptor coverage
#                         rec_pdbs += [pdb_id]
#                     int_pdbs.remove(pdb_id)
#                 if len(int_pdbs) == 0:
#                     break
#             log.info(f"rec_pdbs {rec_pdbs}")

#         if len(rec_pdbs) > 0:
#             # Check ligand coverage
#             try:
#                 lig_allpdb_data = make_request(lig_url, None)
#             except Exception as e:
#                 log.warning(f"Could not make AllPDB request for {ligand_id}, {e}")
#             if lig_allpdb_data:
#                 # Iterate every pdb for the ligand to check if it is one of the receptor-selected ones
#                 for lig_pdb_structure in lig_allpdb_data[ligand_id]["mappings"]:
#                     pdb_id = lig_pdb_structure["entry_id"]
#                     if pdb_id in rec_pdbs:
#                         lig_cov = calculate_coverage(lig_pdb_structure, pdb_id, lig_len)
#                         if lig_cov >= COV_CUTOFF:  # Successful ligand coverage
#                             confirmed_pdbs += [pdb_id]  # It's a valid PDB
#                         rec_pdbs.remove(pdb_id)
#                     if len(rec_pdbs) == 0:
#                         break

#         if len(confirmed_pdbs) == 0:
#             confirmed_pdbs = None
#     return confirmed_pdbs


# def load_negatome_data(negatome_url):
#     """
#     Load the negatome data.

#     Parameters
#     ----------
#     negatome_url : string
#         url to the negatome
#     Returns
#     -------
#     negatome_dic : dict
#         negatome dictionary
#     """
#     negatome_f = negatome_url.split("/")[-1]
#     if not Path(negatome_f).exists():
#         subprocess.run(["wget", negatome_url], check=True)

#     negatome_dic = {}
#     with open(negatome_f, "r") as fh:
#         for line in fh.readlines():
#             partner_a, partner_b = line.split()
#             if partner_a not in negatome_dic.keys():
#                 negatome_dic[partner_a] = []
#             if partner_b not in negatome_dic[partner_a]:
#                 negatome_dic[partner_a] += [partner_b]

#             if partner_b not in negatome_dic.keys():
#                 negatome_dic[partner_b] = []
#             if partner_a not in negatome_dic[partner_b]:
#                 negatome_dic[partner_b] += [partner_a]

#     return negatome_dic


# def check_negatome(uniprot_id, partner_id, negatome_dic):
#     """Checks if the two ids are positive-negative duplicates."""
#     negatome = False
#     if uniprot_id in negatome_dic.keys():
#         if (
#             partner_id in negatome_dic[uniprot_id]
#         ):  # Receptor:Negative partner interaction
#             log.warning(
#                 f"{partner_id} discarded due to positive-negative partner duplicate"
#             )
#             negatome = True
#     if partner_id in negatome_dic.keys():
#         if (
#             uniprot_id in negatome_dic[partner_id]
#         ):  # Negative partner: Receptor interaction
#             log.warning(
#                 f"{partner_id} discarded due to positive-negative partner duplicate"
#             )
#             negatome = True
#     return negatome


# def get_interface(partner, confirmed_int_pdbs):
#     "Get the receptor interface residues for the interaction with a partner"
#     residue_list = {}
#     for resi in partner["residues"]:
#         for int_pdb in resi["interactingPDBEntries"]:
#             int_pdb_id = int_pdb["pdbId"]
#             if int_pdb_id in confirmed_int_pdbs:
#                 if int_pdb_id not in residue_list.keys():
#                     residue_list[int_pdb_id] = []
#                 interface_start = int(resi["startIndex"])
#                 interface_end = int(resi["endIndex"])
#                 if interface_start != interface_end:
#                     # Interface has a few residues.
#                     interface = list(range(interface_start - 1, interface_end))
#                 else:
#                     # Interface is only one residue.
#                     interface = [interface_start]
#                     # Keep track of the residues
#                     for residue in interface:
#                         if residue not in residue_list[int_pdb_id]:
#                             residue_list[int_pdb_id] += [residue]

#     if len(residue_list) == 0:
#         residue_list = None

#     return residue_list


# def check_if_Ig(uniprot_ID, Ig_words):
#     "Use the Uniprot API to find if the given Uniprot ID is an Immunoglobulin."

#     url = f"{PROTEINS_URL}/{uniprot_ID}"
#     Ig = False

#     # Access the API
#     try:
#         protein_data = make_request(url, None)
#     except Exception as e:
#         log.info(f"Could not make InterfaceResidues request for {uniprot_ID}, {e}")
#         return None

#     if protein_data:
#         if len(protein_data) == 0:  # Check if the entry is empty.
#             log.info(f"Proteins API for {uniprot_ID} was empty")
#             return False
#         else:
#             all_names = []  # Store all possible names for the given Uniprot ID
#             if "recommendedName" in protein_data["protein"].keys():  # Main name
#                 try:
#                     main_name = protein_data["protein"]["recommendedName"]["fullName"][
#                         "value"
#                     ]
#                     all_names += [main_name]
#                 except Exception as e:
#                     log.exception(e)
#             if (
#                 "alternativeName" in protein_data["protein"].keys()
#             ):  # Other posible names
#                 try:
#                     for alt_name in protein_data["protein"]["alternativeName"]:
#                         all_names += [alt_name["fullName"]["value"]]
#                 except Exception as e:
#                     log.exception(e)
#             if (
#                 len(all_names) == 0
#             ):  # Something is wrong with this UniprotID API or it's empty
#                 log.info(f"Protein API for {uniprot_ID} did NOT provide a valid name")
#                 return None

#             else:
#                 # Check if any of the Immunoglobulin words is in any of the names that the Uniprot ID has
#                 for Ig_word in Ig_words:
#                     if Ig_word in all_names:
#                         Ig = True
#                         break  # No need to check more
#     return Ig


# def check_multimers_receptor(ligand_id, receptor_id, rec_int_pdbs, outHD):
#     """
#     Check if the receptor protein interacts always as a multimer in a given list of PDBs.

#     Parameters
#     ----------
#     ligand_id : str
#         uniprot ID

#     receptor_id : str
#         uniprot ID

#     rec_int_pdbs : list
#         list of pdb interacting with the receptor

#     outHD : bool
#         if True, exclude heterodimers

#     Returns
#     -------
#     confirmed_pdbs : list
#         list of pdbs that passed the check_multimers assessment
#     """
#     confirmed_pdbs = []  # PDBs to accept
#     # Get info from the API
#     url = f"{INTERFACE_URL}/{ligand_id}"
#     # The request might fail
#     try:
#         interface_data = make_request(url, None)
#     except Exception as e:
#         log.debug(f"Could not make InterfaceResidues request for {ligand_id}, {e}")
#         return False, None

#     if interface_data and len(interface_data) != 0:  # Check if the entry is empty
#         if len(interface_data[ligand_id]) == 0:  # Check if the entry is empty
#             log.debug(f"InterfaceResidues for {ligand_id} was empty")
#         else:
#             # Iterate every interacting partner
#             for partner in interface_data[ligand_id]["data"]:
#                 partner_id = partner["accession"]
#                 if partner_id == ligand_id:  # Ligand is a Homodimer
#                     if outHD:  # Exclude homodimers
#                         log.info(f"{ligand_id} discarded due to being an homodimer")
#                         return False, None

#                 elif partner_id == receptor_id:  # Receptor-partner interaction
#                     for residue in partner["residues"]:  # Iterate every residue
#                         # Check if any of the PDBs for the residue is in the preselected PDBs
#                         for lig_int_pdb in residue["interactingPDBEntries"]:
#                             lig_int_pdb_id = lig_int_pdb["pdbId"]
#                             if (
#                                 lig_int_pdb_id in rec_int_pdbs
#                             ):  # PDB is in the preselected ones
#                                 # Count the number of molecules for the receptor protein in the interacting pdb
#                                 if "chainIds" in lig_int_pdb.keys():
#                                     chain_Ids = lig_int_pdb["chainIds"].split(",")
#                                     if (
#                                         len(chain_Ids) == 1
#                                     ):  # There is only one partner molecule
#                                         if lig_int_pdb_id not in confirmed_pdbs:
#                                             confirmed_pdbs += [
#                                                 lig_int_pdb_id
#                                             ]  # Accept PDB_ID
#                                     rec_int_pdbs.remove(
#                                         lig_int_pdb_id
#                                     )  # Avoid reprocessing same data
#     if len(confirmed_pdbs) == 0:
#         confirmed_pdbs = None
#     return confirmed_pdbs


# def check_only_two(candidate_pdbs):
#     """
#     Check if there are only 2 protein molecules in the PDBs from a given list.

#     Parameters
#     ----------
#     candidate_pdbs : list
#         list of candidate pdb IDs

#     Returns
#     -------
#     confirmed_pdbs : list
#         list of confirmed pdb IDs

#     """
#     confirmed_pdbs = None
#     for pdb in candidate_pdbs:
#         # Get data from the PDB Molecules API
#         pdb_url = f"{PDB_MOLECULES_URL}/{pdb}"
#         try:
#             pdb_data = make_request(pdb_url, None)
#         except Exception as e:
#             log.debug(f"Could not make PDB request for {pdb}, {e}")
#         # if the call is successful
#         if pdb_data:
#             confirmed_pdbs = []
#             n_proteins = 0
#             # Iterate every molecule in the pdb
#             for entity in pdb_data[pdb]:
#                 # Check if entity is a protein, to deal with cofactors and other atoms
#                 if "molecule_type" in entity.keys():
#                     if "polypeptide" in entity["molecule_type"]:  # It is a protein
#                         if "in_struct_asyms" in entity.keys():
#                             n_proteins += len(
#                                 entity["in_struct_asyms"]
#                             )  # Every chain is a different protein
#                             if (
#                                 n_proteins > 2
#                             ):  # Avoid wasting time in processing everything
#                                 break
#             if n_proteins == 2:  # There are only 2 proteins, accept PDB
#                 confirmed_pdbs += [pdb]
#     if len(confirmed_pdbs) == 0:
#         confirmed_pdbs = None
#     return confirmed_pdbs


# def filter_interface_data(uniprot_id, interface_data, filters):
#     """
#     Filters interface data according to the specified filters.

#     Parameters
#     ----------

#     interface_data :

#     filters : dict of booleans
#         employed filters.

#     Returns
#     -------
#     interface_residues : list
#         Interface residues.

#     """
#     receptor_interface_data = {}  # Dictionary to store information
#     # loading negatome data. TODO: generalize this
#     negatome_dic = load_negatome_data(NEGATOME_URL)
#     # place where ProdigyCrystal is going to download the pdbs to proceed with the evaluation
#     # data_folder = None
#     # prodigy_results_pdb = ast.literal_eval(open("/trinity/login/jlopez/interactome/create_benchmark/prodigy_results_pdb.txt", "r").read())
#     # dictionary that stores information of previous ProdigyCrystal runs that tell if
#     # the interactions are biological or not
#     rec_len = len(
#         interface_data[uniprot_id]["sequence"]
#     )  # Lenght of the receptor according to Uniprot
#     log.debug(f"rec_len {rec_len}")
#     # Iterate over every interacting partner
#     accessions = [p["accession"] for p in interface_data[uniprot_id]["data"]]
#     log.debug(f"partners: {accessions}")
#     for partner in interface_data[uniprot_id]["data"]:
#         partner_id = partner["accession"]  # Uniprot_id of the partner
#         log.info(f"Evaluating interaction between {uniprot_id} and {partner_id}")
#         # 1. Check if the interaction is in the Negatome (negative interaction)
#         log.info(
#             f"Checking if {uniprot_id} or {partner_id} are positive-negative partner duplicates"
#         )
#         negatome = check_negatome(
#             uniprot_id, partner_id, negatome_dic
#         )  # true if in negatome
#         if negatome:
#             continue
#         # 2. check if we have to exclude homodimers
#         if filters["outHD"]:  # Exclude Homodimers
#             if partner_id == uniprot_id:  # The receptor itself is a homodimer
#                 log.warning(f"{uniprot_id} discarded due to being a homodimer")
#                 receptor_interface_data = None
#                 # break # Discard the entire receptor. *Jesus: This can be handled in a more thorough manne
#                 continue
#         # 3. check if we have to exclude Immunoglobulins
#         if filters["outIg"]:
#             log.info(f"Checking if {partner_id} is an Immunoglobulin")
#             isIg = check_if_Ig(partner_id, Ig_words)
#             if isIg:  # Interacting partner is an Immunoglobulin
#                 log.warning(f"{partner_id} is an Immunoglobulin, discarding it")
#                 continue
#         # 4. Check if receptor and partner interact as 1 vs 1 in any of the PDB structures
#         if filters["outHM"]:
#             # First step: check if ligand is always a multimer
#             log.info(f"Checking ligand multimers for {partner_id}")
#             excluded_int_pdbs = []  # Multimer pdbs
#             included_int_pdbs = []  # Not multimer pdbs
#             for residue in partner["residues"]:
#                 for int_pdb in residue[
#                     "interactingPDBEntries"
#                 ]:  # Every PDB in which the residue was seen in the interface
#                     int_pdb_id = int_pdb[
#                         "pdbId"
#                     ]  # PDB ID of the structure where the interaction was seen
#                     if (
#                         int_pdb_id not in excluded_int_pdbs
#                         and int_pdb_id not in included_int_pdbs
#                     ):  # Avoid repeated pdbs
#                         # Count the number of molecules for the ligand protein in the interacting pdb
#                         if "chainIds" in int_pdb.keys():
#                             chain_Ids = int_pdb["chainIds"].split(",")
#                             if (
#                                 len(chain_Ids) == 1
#                             ):  # There is only one partner molecule
#                                 included_int_pdbs += [int_pdb_id]  # Accept PDB_ID
#                             else:  # There is more than one partner molecule
#                                 excluded_int_pdbs += [int_pdb_id]  # Discard PDB_ID

#             if (
#                 len(included_int_pdbs) > 0
#             ):  # There are PDBs in which there is only 1 molecule of the partner
#                 # Second step: check if receptor is always a multimer in the included interacting PDBs.
#                 log.info(f"Checking receptor multimers for {uniprot_id}")
#                 confirmed_int_pdbs = check_multimers_receptor(
#                     partner_id, uniprot_id, included_int_pdbs, filters["outHD"]
#                 )
#                 if confirmed_int_pdbs is None:  # There are 1:1 PPI interactions
#                     # (Optional step): check if only two protein molecules in the pdb structure
#                     if filters["t2m"]:
#                         log.info("Checking only 2 proteins in the pdb files.")
#                         confirmed_int_pdbs = check_only_two(confirmed_int_pdbs)
#                         if confirmed_int_pdbs is None:
#                             log.warning(
#                                 f"{uniprot_id}_{partner_id} discarded due to no pdb with only 2 proteins"
#                             )
#                             continue
#                 else:
#                     log.warning(
#                         f"{uniprot_id}_{partner_id} discarded due to receptor only interacts as multimer"
#                     )
#                     continue
#             else:
#                 log.warning(
#                     f"{uniprot_id}_{partner_id} discarded due to ligand only interacts as multimer"
#                 )
#                 continue
#         else:
#             confirmed_int_pdbs = partner["allPDBEntries"]
#         # 5. check if coverage of the proteins in the interacting PDBs is over the established threshold
#         if filters["cpc"]:
#             log.info(f"Checking partners coverage for {uniprot_id}_{partner_id}")
#             confirmed_int_pdbs = check_partners_coverage(
#                 uniprot_id, rec_len, partner_id, confirmed_int_pdbs
#             )
#             if confirmed_int_pdbs is None:
#                 log.warning(
#                     f"{uniprot_id}_{partner_id} discarded due to coverage under cutoff"
#                 )
#                 continue
#         # 6. check: if the interactions for the confirmed pdbs are biological by Prodigy Crystal
#         # if filters["cbi"]: # Check: if the interactions for the confirmed pdbs are biological by Prodigy Crystal
#         #     log.info(f"Checking biological interaction for {uniprot_id}-{partner_id}")
#         #     confirmed_int_pdbs, prodigy_results_pdb = validate_interaction(uniprot_id, partner_id, confirmed_int_pdbs, prodigy_results_pdb, data_folder)
#         #     if confirmed_int_pdbs == None:
#         #         log.warning(f"{partner_id} discarded due to no biological interaction")
#         #         continue
#         if filters["cbi"]:
#             log.warning("Checking biological interaction not implemented yet.")
#         # Get interface residues of the single PPI interaction
#         log.info(f"Checking matching interfaces in the structure of {uniprot_id}")
#         interface = get_interface(partner, confirmed_int_pdbs)
#         if interface is not None:  # There are valid PDBs and valid residue
#             # Check the interface of each PDB separatedly
#             for int_pdb in interface.keys():
#                 inter_check = True
#                 # Check if interface residues are present in the receptor structure
#                 # for residue in interface[int_pdb]: # now rec_residues is not used anymore
#                 #    if residue not in rec_residues:
#                 #        inter_check = False
#                 #        break
#                 if inter_check:  # The interface is valid
#                     if partner_id not in receptor_interface_data:
#                         log.info("It's a match")
#                         receptor_interface_data[partner_id] = []
#                     # Add residues to partner interface
#                     for residue in interface[int_pdb]:
#                         if residue not in receptor_interface_data[partner_id]:
#                             receptor_interface_data[partner_id] += [residue]

#             if (
#                 partner_id not in receptor_interface_data.keys()
#             ):  # No valid interfaces in receptor structure
#                 # log.warning(f"{partner_id} discarded due to interface with {receptor_id}_{receptor_pdb} not in receptor structure")
#                 log.warning(
#                     f"{partner_id} discarded due to interface with {uniprot_id} not valid"
#                 )
#         else:
#             # log.warning(f"{partner_id} discarded due to interface with {receptor_id}_{receptor_pdb} not valid")
#             log.warning(
#                 f"{partner_id} discarded due to interface with {uniprot_id} not valid"
#             )

#     if (
#         receptor_interface_data is not None
#     ):  # In outHD it can be None if receptor is a homodimer
#         if len(receptor_interface_data) == 0:
#             receptor_interface_data = None

#     # return receptor_interface_data, prodigy_results_pdb
#     return receptor_interface_data
