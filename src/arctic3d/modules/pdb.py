import logging
from pathlib import Path

from pdbtools.pdb_fetch import fetch_structure
from pdbtools.pdb_selaltloc import select_by_occupancy
from pdbtools.pdb_selchain import select_chain
from pdbtools.pdb_tidy import tidy_pdbfile

from arctic3d.functions import make_request

log = logging.getLogger("arctic3dlog")

BESTPDB_URL = "https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures"


def fetch_pdb(pdb_id):
    """
    Fetch PDB file.
    """
    log.debug(f"Fetching PDB file {pdb_id}")
    out_pdb_fname = Path(f"{pdb_id}.pdb")
    with open(out_pdb_fname, "w") as f:
        for line in fetch_structure(pdb_id):
            f.write(line)
    return out_pdb_fname


def selchain_pdb(inp_pdb_f, chain):
    """
    Select chain from PDB file.
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
    """
    log.debug("Removing non-ATOM lines from PDB file")
    out_pdb_fname = Path(f"{inp_pdb_f.stem}-atoms.pdb")
    with open(inp_pdb_f, "r") as pdb_fh:
        with open(out_pdb_fname, "w") as f:
            for line in pdb_fh:
                if line.startswith("ATOM"):
                    f.write(line)
    return out_pdb_fname


def get_best_pdb(uniprot_id):
    """
    Get best PDB ID.
    """
    pdb_dict = {}
    url = f"{BESTPDB_URL}/{uniprot_id}"
    try:
        pdb_dict = make_request(url, None)
    except Exception as e:
        log.warning(f"Could not make BestStructure request for {uniprot_id}, {e}")
        return pdb_dict

    # the first entry will be the one with highest coverage/lower resolution
    top_hit = pdb_dict[uniprot_id][0]

    pdb_id = top_hit["pdb_id"]
    chain_id = top_hit["chain_id"]
    # TODO: Add a check for minimum coverage/resolution
    coverage = top_hit["coverage"]
    resolution = top_hit["resolution"]
    start = top_hit["unp_start"]
    end = top_hit["unp_end"]

    log.info(
        f"BestPDB hit for {uniprot_id}: {pdb_id}_{chain_id} {coverage:.2f} coverage {resolution:.2f} Angstrom / start {start} end {end}"
    )

    pdb_f = fetch_pdb(pdb_id)
    atoms_pdb_f = keep_atoms(pdb_f)
    chained_pdb_f = selchain_pdb(atoms_pdb_f, chain_id)
    occ_pdb_f = occ_pdb(chained_pdb_f)
    tidy_pdb_f = tidy_pdb(occ_pdb_f)

    tidy_pdb_f.rename(f"{uniprot_id}-{pdb_id}-{chain_id}.pdb")

    pdb_f.unlink()
    atoms_pdb_f.unlink()
    chained_pdb_f.unlink()
    occ_pdb_f.unlink()

    return tidy_pdb_f
