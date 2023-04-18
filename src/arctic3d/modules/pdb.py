import logging
import os
from pathlib import Path

import jsonpickle
import MDAnalysis as mda
import requests
from pdbecif.mmcif_io import MMCIF2Dict
from pdbtools.pdb_selaltloc import select_by_occupancy
from pdbtools.pdb_selchain import select_chain
from pdbtools.pdb_tidy import tidy_pdbfile
from pdbtools.pdb_selmodel import select_model


from arctic3d.functions import make_request
from arctic3d.modules.interface_matrix import filter_interfaces

log = logging.getLogger("arctic3d.log")

BESTPDB_URL = "https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures"
PDBE_URL = "https://www.ebi.ac.uk/pdbe/entry-files/download"


def fetch_updated_cif(pdb_id, cif_fname):
    """
    Fetch updated cif from PDBE database.

    Parameters
    ----------
    pdb_id : str
        PDB ID
    cif_fname : str or Path
        name of the output cif file
    """
    log.debug(f"Fetching updated CIF file {pdb_id} from PDBE")
    cif_response = requests.get(f"{PDBE_URL}/{pdb_id}_updated.cif")
    if cif_response.status_code != 200:
        log.warning(f"Could not fetch CIF file for {pdb_id}")
        return None

    with open(cif_fname, "wb") as wfile:
        wfile.write(cif_response.content)


def get_cif_dict(cif_name):
    """
    Convert cif file to dict.

    Parameters
    ----------
    cif_name : str or Path
        cif filename

    Returns
    -------
    cif_dict : dict
        cif dictionary
    """
    mmcif_dict = MMCIF2Dict()
    cif_dict = mmcif_dict.parse(cif_name)
    return cif_dict


def get_numbering_dict(pdb_id, cif_dict, uniprot_id, chain_id, key="pdb"):
    """
    gets the numbering correspondence between the pdb file and the uniprot
    sequence from the cif dict.

    Parameters
    ----------
    pdb_id : str
        PDB ID
    cif_dict : dict
        cif dictionary
    uniprot_id : str
        uniprot ID to be used (many IDs may exist in the .cif file)
    chain_id : str
        chain ID to be used
    key : str
        key to use for the numbering dict, either "uniprot" or "pdb"

    Returns
    -------
    numbering_dict : dict
        pdb-resid : key-value dictionary
            Example (key=pdb) : {"GLY-A-16" : 20, "TYR-A-17" : 21, ... }
            Example (key=uniprot) : {20 : "GLY-A-16", 21 : "TYR-A-17", ... }
    """
    atomsite_dict = cif_dict[pdb_id.upper()]["_atom_site"]
    numbering_dict = {}
    prev_residue_key = None
    len_sifts_mapping = len(atomsite_dict["auth_seq_id"])
    for resid in range(len_sifts_mapping):
        if (
            atomsite_dict["pdbx_sifts_xref_db_acc"][resid] == uniprot_id
            and atomsite_dict["auth_asym_id"][resid] == chain_id
        ):
            # check for an insertion code
            ins_code = atomsite_dict["pdbx_PDB_ins_code"][resid]
            residue_key = (
                f"{atomsite_dict['auth_comp_id'][resid]}"
                f"-{atomsite_dict['auth_asym_id'][resid]}"
                f"-{atomsite_dict['auth_seq_id'][resid]}"
                f"{('' if ins_code == '?' else f'-{ins_code}')}"
            )
            unp_num = atomsite_dict["pdbx_sifts_xref_db_num"][resid]
            if residue_key != prev_residue_key:  # not a duplicate entry
                if key == "pdb":
                    numbering_dict[residue_key] = unp_num
                elif key == "uniprot":
                    numbering_dict[unp_num] = residue_key
                else:
                    raise ValueError(f"key {key} not recognized")
                prev_residue_key = residue_key
    # log.debug(f"numbering dict {numbering_dict}")
    return numbering_dict


def renumber_interfaces_from_cif(
    pdb_id, uniprot_id, chain_id, interface_residues
):
    """
    Renumbers a list of interfaces based on the information coming from the
    corresponding updated cif file.

    Parameters
    ----------
    pdb_id : str
        PDB ID
    uniprot_id : str
        uniprot ID to be used
    chain_id : str
        chain ID to be used
    interfaces_residues : list
        list of interfaces residues
    """

    cif_fname = Path(f"{pdb_id}_updated.cif")
    if not cif_fname.is_file():
        fetch_updated_cif(pdb_id, cif_fname)
    cif_dict = get_cif_dict(cif_fname)

    # retrieve mapping
    numbering_dict = get_numbering_dict(
        pdb_id, cif_dict, uniprot_id, chain_id, key="uniprot"
    )
    # log.debug(f"numbering_dict {numbering_dict}")
    if any(numbering_dict):
        unique_resids = set(
            value for values in interface_residues.values() for value in values
        )
        renum_residues = {}  # dictionary of renumbered residues
        for residue in unique_resids:
            str_res = str(residue)
            if str_res in numbering_dict.keys():
                # log.debug(f"Residue {residue} not found in cif file")
                int_residue = int(numbering_dict[str_res].split("-")[2])
                renum_residues[residue] = int_residue
            else:
                # log.debug(f"Residue {residue} not found in cif file")
                renum_residues[residue] = None
        # renumbering interfaces
        renum_interfaces = {}
        for interface, residues in interface_residues.items():
            renum_residues_list = []
            for residue in residues:
                if residue is not None:
                    renum_residues_list.append(renum_residues[residue])
            renum_interfaces[interface] = renum_residues_list
    else:
        log.info(f"Renumbering failed for pdb {pdb_id}-{chain_id}")
        renum_interfaces = None
    # log.debug(f"renum_interfaces {renum_interfaces}")
    return renum_interfaces, cif_fname


def renumber_pdb_from_cif(pdb_id, uniprot_id, chain_id, pdb_fname):
    """
    Renumbers a pdb file based on the information coming from the corresponding
    updated cif file.

    Parameters
    ----------
    pdb_id : str
        PDB ID
    chain_id : str
        chain ID to be used
    uniprot_id : str
        uniprot ID to be used
    pdb_fname : str or Path
        input pdb file

    Returns
    -------
    pdb_renum_fname : Path
        renumbered pdb filename
    """
    cif_fname = Path(f"{pdb_id}_updated.cif")
    if not cif_fname.is_file():
        fetch_updated_cif(pdb_id, cif_fname)
    cif_dict = get_cif_dict(cif_fname)

    # retrieve mapping
    numbering_dict = get_numbering_dict(
        pdb_id, cif_dict, uniprot_id, chain_id, key="pdb"
    )

    # we do not check if all residues in pdb_fname have
    #   been correctly renumbered
    # we only check it's not empty (it could be empty
    #   if the cif does not contain
    # the uniprot information)
    if any(numbering_dict):
        log.info(f"Renumbering pdb {pdb_fname}")
        pdb_renum_fname = Path(f"{pdb_fname.stem}_renum.pdb")

        records = ("ATOM", "TER")
        file_content = ""
        with open(pdb_renum_fname, "w") as wfile:
            with open(pdb_fname, "r") as rfile:
                for ln in rfile:
                    if ln.startswith(records):
                        resid = ln[22:26].strip()
                        ins_code = "-" + ln[26] if ln[26] != " " else ""
                        residue_key = (
                            f"{ln[17:20].strip()}"
                            f"-{ln[20:22].strip()}"
                            f"-{resid}{ins_code}"
                        )  # resname-chain_id-resid-ins_code

                        # the residues in the pdb_fname that do not have an
                        #   entry in the numbering_dict
                        # are discarded. It may happen that the same chain
                        #   in the input pdb is associated to several
                        # uniprot ids (especially in old files)
                        if residue_key in numbering_dict.keys():
                            n_spaces = 4 - len(
                                str(numbering_dict[residue_key])
                            )
                            # there's always one space after to remove
                            #   alternate occupancies
                            resid_str = (
                                f"{' ' * n_spaces}"
                                f"{numbering_dict[residue_key]} "
                            )
                            file_content += f"{ln[:22]}{resid_str}{ln[27:]}"
                    else:
                        file_content += f"{ln}"
                wfile.write(file_content)
    else:
        log.info(f"Renumbering failed for pdb {pdb_fname}")
        pdb_renum_fname = None
    return pdb_renum_fname, cif_fname


def fetch_pdb_files(pdb_to_fetch):
    """
    Fetches the pdb files from PDBe database.

    Parameters
    ----------
    pdb_to_fetch : list
        list of pdb hits to fetch

    Returns
    -------
    validated_pdbs : list
        list of tuples (pdb_file, hit)
    """
    validated_pdbs = []
    valid_pdb_set = set()  # set of valid pdb IDs
    for hit in pdb_to_fetch:
        pdb_id = hit["pdb_id"]
        pdb_fname = f"{pdb_id}.pdb"
        if pdb_fname not in os.listdir():
            pdb_f = fetch_pdb(pdb_id)
        else:
            pdb_f = Path(pdb_fname)
        if pdb_f is not None:
            validated_pdbs.append((pdb_f, hit))
            if pdb_id not in valid_pdb_set:
                valid_pdb_set.add(pdb_id)
    return validated_pdbs


def fetch_pdb(pdb_id):
    """
    Fetches the pdb from PDBe database.

    This is the un-renumbered pdb.

    Parameters
    ----------
    pdb_id : str
        pdb target id
    Returns
    -------
    out_pdb_fname : Path
        pdb filename
    """
    log.debug(f"Fetching PDB file {pdb_id} from PDBE")
    response = requests.get(f"{PDBE_URL}/pdb{pdb_id}.ent")
    if response.status_code != 200:
        log.warning(f"Could not fetch PDB file for {pdb_id}")
        return None

    out_pdb_fname = Path(f"{pdb_id}.pdb")
    with open(out_pdb_fname, "wb") as wfile:
        wfile.write(response.content)

    return out_pdb_fname


def selchain_pdb(inp_pdb_f, chain):
    """
    Select chain from PDB file.

    Parameters
    ----------
    inp_pdb_f : Path
        Path to PDB file.
    chain : str
        Chain ID.

    Returns
    -------
    out_pdb_fname : Path
        Path to PDB file.
    """
    # log.debug(f"Selecting chain {chain} from PDB file")
    out_pdb_fname = Path(f"{inp_pdb_f.stem}-{chain}.pdb")
    with open(inp_pdb_f, "r") as pdb_fh:
        with open(out_pdb_fname, "w") as f:
            for line in select_chain(pdb_fh, chain):
                f.write(line)
    return out_pdb_fname


def selmodel_pdb(inp_pdb_f, model_id=1):
    """
    Select model from PDB file.

    Parameters
    ----------
    inp_pdb_f : Path
        Path to PDB file.
    model_id : int, optional
        Model ID, by default 1

    Returns
    -------
    out_pdb_fname : Path
        Path to PDB file.
    """
    # log.debug(f"Selecting model {model_id} from PDB file")
    out_pdb_fname = Path(f"{inp_pdb_f.stem}-model{model_id}.pdb")
    with open(inp_pdb_f, "r") as pdb_fh:
        with open(out_pdb_fname, "w") as f:
            line = ""
            for line in select_model(pdb_fh, [model_id]):
                f.write(line)
    return out_pdb_fname


def tidy_pdb(inp_pdb_f):
    """
    Tidy PDB file.

    Parameters
    ----------
    inp_pdb_f : Path
        Path to PDB file.

    Returns
    -------
    Path
        Path to PDB file.
    """
    # log.debug("Tidying PDB file")
    out_pdb_fname = Path(f"{inp_pdb_f.stem}-tidy.pdb")
    with open(inp_pdb_f, "r") as pdb_fh:
        with open(out_pdb_fname, "w") as f:
            for line in tidy_pdbfile(pdb_fh):
                f.write(line)
    return out_pdb_fname


def occ_pdb(inp_pdb_f):
    """
    Select residues with highest occupancy.

    Parameters
    ----------
    inp_pdb_f : Path
        Path to PDB file.

    Returns
    -------
    out_pdb_fname : Path
        Path to PDB file.
    """
    # log.debug("Selecting residues with highest occupancy")
    out_pdb_fname = Path(f"{inp_pdb_f.stem}-occ.pdb")
    with open(inp_pdb_f, "r") as pdb_fh:
        with open(out_pdb_fname, "w") as f:
            for line in select_by_occupancy(pdb_fh):
                f.write(line)
    return out_pdb_fname


def keep_atoms(inp_pdb_f):
    """
    Keep atoms.

    Parameters
    ----------
    inp_pdb_f : Path
        Path to PDB file.

    Returns
    -------
    out_pdb_fname : Path
        Path to PDB file.
    """
    # log.debug(f"Removing non-ATOM lines from PDB file {inp_pdb_f}")
    out_pdb_fname = Path(f"{inp_pdb_f.stem}-atoms.pdb")
    with open(inp_pdb_f, "r") as pdb_fh:
        with open(out_pdb_fname, "w") as f:
            for line in pdb_fh:
                if line.startswith("ATOM"):
                    f.write(line)
    return out_pdb_fname


def validate_api_hit(
    fetch_list,
    resolution_cutoff=4.0,
    coverage_cutoff=0.0,
    max_pdb_num=20,
):
    """
    Validate PDB fetch request file.

    Parameters
    ----------
    fetch_list : list
        List containing dictionaries of hits.
    pdb_renum_db : str or Path or None
        path to the pdb renum local db
    resolution_cutoff : float
        Resolution cutoff.
    coverage_cutoff : float
        Coverage cutoff.
    max_pdb_num : int
        Maximum number of pdb files to fetch.

    Returns
    -------
    validated_pdbs : list
        List of (pdb_f, hit) tuples
    """
    pdbs_to_fetch = []
    for hit in fetch_list:
        check_list = []
        pdb_id = hit["pdb_id"]
        coverage = hit["coverage"]
        resolution = hit["resolution"]
        exp_method = hit["experimental_method"]

        # check coverage value
        if coverage > coverage_cutoff:
            check_list.append(True)
        else:
            check_list.append(False)
        # check resolution value
        if resolution is None:
            if "NMR" in exp_method:
                check_list.append(True)
            else:
                check_list.append(False)
        elif resolution < resolution_cutoff:
            check_list.append(True)
        else:
            check_list.append(False)

        if all(check_list):
            pdbs_to_fetch.append(hit)
        else:
            log.debug(f"{pdb_id} failed validation")
    log.info(f"Found {len(pdbs_to_fetch)} valid PDBs to fetch")
    # downloading a list of good pdbs
    validated_pdbs = fetch_pdb_files(pdbs_to_fetch[:max_pdb_num])
    log.info(f"Fetched {len(validated_pdbs)} valid PDBs")
    return validated_pdbs


def preprocess_pdb(pdb_fname, chain_id):
    """
    Apply a set of transformations to an input pdb file.

    Parameters
    ----------
    pdb_fname : str or Path
        input pdb file
    chain_id : str
        chain ID

    Returns
    -------
    tidy_pdb_f : Path
        preprocessed pdb file
    """
    model_pdb_f = selmodel_pdb(pdb_fname)
    atoms_pdb_f = keep_atoms(model_pdb_f)
    chained_pdb_f = selchain_pdb(atoms_pdb_f, chain_id)
    occ_pdb_f = occ_pdb(chained_pdb_f)
    tidy_pdb_f = tidy_pdb(occ_pdb_f)

    atoms_pdb_f.unlink()
    chained_pdb_f.unlink()
    occ_pdb_f.unlink()

    return tidy_pdb_f


def unlink_files(suffix="pdb", to_exclude=None):
    """
    Remove all files with suffix in the cwd except for those in to_exclude.

    Parameters
    ----------
    suffix : str

    to_exclude : None or list
        files to exclude
    """
    suffix_fnames = list(Path(".").glob(f"*{suffix}"))
    for fname in suffix_fnames:
        fpath = Path(fname)
        if fpath.is_file() and fpath not in to_exclude:
            fpath.unlink()


def get_maxint_pdb(
    validated_pdbs, interface_residues, uniprot_id, numbering="pdb"
):
    """
    Get PDB ID that retains the most interfaces.

    Parameters
    ----------
    validated_pdbs : list
        List of (pdb_f, hit) tuples
    interface_residues : dict
        Dictionary of all the interfaces (each one with its uniprot ID as key)
    uniprot_id : str
        Uniprot ID
    numbering : str
        what to renumber? 'pdb' for pdb files, 'resi' for interface residues
    """
    log.info("Selecting pdb retaining the most interfaces")
    cif_f, pdb_f, hit, filtered_interfaces = None, None, None, None
    if validated_pdbs != []:
        max_nint = 0
        for curr_pdb, curr_hit in validated_pdbs:
            chain_id = curr_hit["chain_id"]
            pdb_id = curr_hit["pdb_id"]

            # refactor renumbering
            tidy_pdb_f = preprocess_pdb(curr_pdb, chain_id)

            if numbering == "pdb":  # renumber the pdb files
                curr_pdb_f, curr_cif_f = renumber_pdb_from_cif(
                    pdb_id, uniprot_id, chain_id, tidy_pdb_f
                )
                curr_interface_residues = interface_residues
            elif numbering == "resi":  # renumber the interface residues
                curr_pdb_f = tidy_pdb_f
                (
                    curr_interface_residues,
                    curr_cif_f,
                ) = renumber_interfaces_from_cif(
                    pdb_id, uniprot_id, chain_id, interface_residues
                )
            else:
                raise ValueError(f"Unknown numbering option: {numbering}")
            # load pdb file. If there is an error, skip to the next one
            try:
                mdu = mda.Universe(curr_pdb_f)
            except Exception as e:
                log.error(f"Error loading {curr_pdb_f}: {e}")
                continue

            selection_string = f"name CA and chainID {chain_id}"
            pdb_resids = mdu.select_atoms(selection_string).resids
            tmp_filtered_interfaces = filter_interfaces(
                curr_interface_residues, pdb_resids
            )
            curr_nint = len(tmp_filtered_interfaces)
            if curr_nint > max_nint:  # update "best" hit
                max_nint = curr_nint
                filtered_interfaces = tmp_filtered_interfaces.copy()
                pdb_f = curr_pdb_f
                cif_f = curr_cif_f
                hit = curr_hit
        # unlink pdb files
        unlink_files("pdb", to_exclude=[pdb_f])
        unlink_files("cif", to_exclude=[cif_f])

        if max_nint != 0:
            log.info(f"filtered_interfaces {filtered_interfaces}")
            log.info(f"pdb {pdb_f} retains the most interfaces ({max_nint})")
    return pdb_f, cif_f, hit, filtered_interfaces


def filter_pdb_list(fetch_list, pdb_to_use=None, chain_to_use=None):
    """
    Filter the PDB fetch list.

    Parameters
    ----------
    fetch_list : list
        List containing dictionaries of hits.
    pdb_to_use : str
        Pdb code to be used.

    Returns
    -------
    reduced_list : list
        List containing only the pdb_to_use hit
    """

    reduced_list = []
    for hit in fetch_list:
        pdb_id = hit["pdb_id"]
        chain_id = hit["chain_id"]
        pdb_check, chain_check = True, True

        if pdb_to_use and pdb_id != pdb_to_use:
            pdb_check = False
        if chain_to_use and chain_id != chain_to_use:
            chain_check = False
        if (pdb_check, chain_check) == (True, True):
            reduced_list.append(hit)

    if len(reduced_list) == 0:
        log.warning(f"PDB ID {pdb_to_use} not found in fetched pdb list.")
    return reduced_list


def get_best_pdb(
    uniprot_id,
    interface_residues,
    pdb_to_use=None,
    chain_to_use=None,
    pdb_data=None,
    numbering="pdb",
):
    """
    Get best PDB ID.

    Parameters
    ----------
    uniprot_id : str
        Uniprot ID.
    interface_residues : dict
        Dictionary of all the interfaces (each one with its uniprot ID as key).
    pdb_to_use : str (default None)
        Pdb code to be used.
    chain_to_use : str (default None)
        Chain id to be used.
    pdb_data : Path or None
        pdb json file for offline mode.
    numbering : str (default pdb)
        what to renumber, either the pdb files or the interface residues

    Returns
    -------
    Path or None
        Path to PDB file or None if no PDB file was found.
    filtered_interfaces : dict or None
        Dictionary of the retained and filtered interfaces.
    """
    pdb_dict = {}
    if not pdb_data:
        url = f"{BESTPDB_URL}/{uniprot_id}"
        try:
            pdb_dict = make_request(url, None)
        except Exception as e:
            log.warning(
                f"Could not make BestStructure request for {uniprot_id}, {e}"
            )
            return
    else:
        try:
            pdb_dict = jsonpickle.decode(open(pdb_data, "r").read())
        except Exception as e:
            log.warning(f"Could not read input interface_data {pdb_data}, {e}")
            return

    # if pdb_to_use is not None, already filter the list
    if pdb_to_use:
        pdb_to_use = pdb_to_use.lower()
    if chain_to_use:
        chain_to_use = chain_to_use.upper()
    pdb_list = filter_pdb_list(pdb_dict[uniprot_id], pdb_to_use, chain_to_use)

    validated_pdbs = validate_api_hit(pdb_list)

    pdb_f, cif_f, top_hit, filtered_interfaces = get_maxint_pdb(
        validated_pdbs, interface_residues, uniprot_id, numbering=numbering
    )

    if pdb_f is None or cif_f is None:
        log.warning(f"Could not fetch PDB/mmcif file for {uniprot_id}")
        return None, None

    pdb_id = top_hit["pdb_id"]
    chain_id = top_hit["chain_id"]
    coverage = top_hit["coverage"]
    resolution = top_hit["resolution"]
    start = top_hit["unp_start"]
    end = top_hit["unp_end"]

    log.info(
        f"BestPDB hit for {uniprot_id}:"
        f" {pdb_id}_{chain_id} {coverage} coverage"
        f" {resolution} Angstrom / start {start} end {end}"
    )

    processed_pdb = pdb_f.rename(f"{uniprot_id}-{pdb_id}-{chain_id}.pdb")

    return processed_pdb, cif_f, filtered_interfaces
