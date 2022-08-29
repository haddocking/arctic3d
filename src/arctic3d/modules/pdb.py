import gzip
import logging
import os
import tempfile
from pathlib import Path

import MDAnalysis as mda
import requests
from pdbtools.pdb_selaltloc import select_by_occupancy
from pdbtools.pdb_selchain import select_chain
from pdbtools.pdb_tidy import tidy_pdbfile

from arctic3d.functions import make_request
from arctic3d.modules.interface_matrix import filter_interfaces

log = logging.getLogger("arctic3dlog")

BESTPDB_URL = "https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures"
PDBRENUM_URL = "http://dunbrack3.fccc.edu/PDBrenum/output_PDB"


def fetch_pdbrenum(pdb_id):
    """
    Fetch PDB file.

    Parameters
    ----------
    pdb_id : str
        PDB ID.

    Returns
    -------
    out_pdb_fname : Path
        Path to PDB file.
    """
    log.debug(f"Fetching PDB file {pdb_id} from PDBrenum")
    response = requests.get(f"{PDBRENUM_URL}/{pdb_id}_renum.pdb.gz")
    if response.status_code != 200:
        log.warning(f"Could not fetch PDB file for {pdb_id}")
        return None

    temp_gz = tempfile.NamedTemporaryFile(mode="wb", delete=False, suffix=".gz")
    temp_gz.write(response.content)
    temp_gz.close()

    out_pdb_fname = Path(f"{pdb_id}.pdb")

    with open(out_pdb_fname, "w") as f:
        with gzip.open(temp_gz.name, "rt") as gz:
            for line in gz:
                f.write(line)

    Path(temp_gz.name).unlink()

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
    log.debug(f"Selecting chain {chain} from PDB file")
    out_pdb_fname = Path(f"{inp_pdb_f.stem}-{chain}.pdb")
    with open(inp_pdb_f, "r") as pdb_fh:
        with open(out_pdb_fname, "w") as f:
            for line in select_chain(pdb_fh, chain):
                f.write(line)
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
    log.debug("Tidying PDB file")
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
    log.debug("Selecting residues with highest occupancy")
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
    log.debug("Removing non-ATOM lines from PDB file")
    out_pdb_fname = Path(f"{inp_pdb_f.stem}-atoms.pdb")
    with open(inp_pdb_f, "r") as pdb_fh:
        with open(out_pdb_fname, "w") as f:
            for line in pdb_fh:
                if line.startswith("ATOM"):
                    f.write(line)
    return out_pdb_fname


def validate_api_hit(
    fetch_list, resolution_cutoff=3.0, coverage_cutoff=0.7, max_pdb_renum=20
):
    """
    Validate PDB fetch request file.

    Parameters
    ----------
    fetch_dic : list
        List containing dictionaries of hits.
    resolution_cutoff : float
        Resolution cutoff.
    coverage_cutoff : float
        Coverage cutoff.
    max_pdb_renum : int
        Maximum number of pdb files to fetch.

    Returns
    -------
    validated_pdbs : list
        List of (pdb_f, hit) tuples
    """
    validated_pdbs = []  # list of good pdbs
    valid_pdb_set = set()  # set of valid pdb IDs

    for hit in fetch_list[:max_pdb_renum]:
        check_list = []
        pdb_id = hit["pdb_id"]
        coverage = hit["coverage"]
        resolution = hit["resolution"]

        pdb_f = fetch_pdbrenum(pdb_id)
        if pdb_f is not None:
            check_list.append(True)
        else:
            check_list.append(False)

        if coverage > coverage_cutoff:
            check_list.append(True)
        else:
            check_list.append(False)

        if resolution is None:
            check_list.append(False)
        elif resolution < resolution_cutoff:
            check_list.append(True)
        else:
            check_list.append(False)

        if all(check_list):
            validated_pdbs.append((pdb_f, hit))
            if pdb_id not in valid_pdb_set:
                valid_pdb_set.add(pdb_id)
        else:
            log.debug(f"{pdb_id} failed validation")
            if (
                pdb_f is not None and pdb_id not in valid_pdb_set
            ):  # pdb_f could be None or the pdb (another chain) could be valid and the file should not be removed
                os.unlink(pdb_f)
    return validated_pdbs


def get_maxint_pdb(validated_pdbs, interface_residues):
    """
    Get PDB ID that retains the most interfaces.

    Parameters
    ----------
    validated_pdbs : list
        List of (pdb_f, hit) tuples
    interface_residues : dict
        Dictionary of all the interfaces (each one with its uniprot ID as key)

    Returns
    -------
    pdb_f : Path or None
        Path to PDB file.
    hit : dict or None
        Interface API hit.
    filtered_interfaces : dict or None
        Dictionary of the retained and filtered interfaces.
    """
    pdb_f, hit, filtered_interfaces = None, None, None
    if validated_pdbs != []:
        max_nint = 0
        for curr_pdb, curr_hit in validated_pdbs:
            mdu = mda.Universe(curr_pdb)
            pdb_resids = mdu.select_atoms("name CA").resids
            tmp_filtered_interfaces = filter_interfaces(interface_residues, pdb_resids)
            curr_nint = len(tmp_filtered_interfaces)
            if curr_nint > max_nint:  # update "best" hit
                max_nint = curr_nint
                filtered_interfaces = tmp_filtered_interfaces.copy()
                pdb_f = curr_pdb
                hit = curr_hit
        # unlink pdb files
        for curr_pdb, curr_hit in validated_pdbs:
            if os.path.exists(curr_pdb):
                if curr_pdb != pdb_f:
                    os.unlink(curr_pdb)
        if max_nint != 0:
            log.info(f"filtered_interfaces {filtered_interfaces}")
            log.info(f"pdb {pdb_f} retains the most interfaces ({max_nint})")
    return pdb_f, hit, filtered_interfaces


def get_best_pdb(uniprot_id, interface_residues):
    """
    Get best PDB ID.

    Parameters
    ----------
    uniprot_id : str
        Uniprot ID.
    interface_residues : dict
        Dictionary of all the interfaces (each one with its uniprot ID as key).

    Returns
    -------
    Path or None
        Path to PDB file or None if no PDB file was found.
    filtered_interfaces : dict or None
        Dictionary of the retained and filtered interfaces.
    """
    pdb_dict = {}
    url = f"{BESTPDB_URL}/{uniprot_id}"
    try:
        pdb_dict = make_request(url, None)
    except Exception as e:
        log.warning(f"Could not make BestStructure request for {uniprot_id}, {e}")
        return

    validated_pdbs = validate_api_hit(pdb_dict[uniprot_id])

    pdb_f, top_hit, filtered_interfaces = get_maxint_pdb(
        validated_pdbs, interface_residues
    )

    if pdb_f is None:
        log.warning(f"Could not fetch PDB file for {uniprot_id}")
        return None, None

    pdb_id = top_hit["pdb_id"]
    chain_id = top_hit["chain_id"]
    coverage = top_hit["coverage"]
    resolution = top_hit["resolution"]
    start = top_hit["unp_start"]
    end = top_hit["unp_end"]

    log.info(
        f"BestPDB hit for {uniprot_id}: {pdb_id}_{chain_id} {coverage:.2f} coverage {resolution:.2f} Angstrom / start {start} end {end}"
    )

    atoms_pdb_f = keep_atoms(pdb_f)
    chained_pdb_f = selchain_pdb(atoms_pdb_f, chain_id)
    occ_pdb_f = occ_pdb(chained_pdb_f)
    tidy_pdb_f = tidy_pdb(occ_pdb_f)

    pdb_f.unlink()
    atoms_pdb_f.unlink()
    chained_pdb_f.unlink()
    occ_pdb_f.unlink()

    processed_pdb = tidy_pdb_f.rename(f"{uniprot_id}-{pdb_id}-{chain_id}.pdb")

    return processed_pdb, filtered_interfaces
