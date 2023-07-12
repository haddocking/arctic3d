import logging
import os
from pathlib import Path

import jsonpickle
import MDAnalysis as mda
import requests
from pdbecif.mmcif_io import MMCIF2Dict

# from pdbtools.pdb_selaltloc import select_by_occupancy
from pdbtools.pdb_selchain import select_chain
from pdbtools.pdb_tidy import tidy_pdbfile
from pdbtools.pdb_selmodel import select_model


from arctic3d.functions import make_request
from arctic3d.modules.interface_matrix import filter_interfaces

log = logging.getLogger("arctic3d.log")

BESTPDB_URL = "https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures"
PDBE_URL = "https://www.ebi.ac.uk/pdbe/entry-files/download"


def _remove_altloc(lines):
    # the altloc ID is removed in processed altloc lines
    for line_num, line in lines:
        yield (line_num, line[:16] + " " + line[17:])


def _flush(register, option, others):
    """
    Processes the collected atoms according to the selaltloc option.
    """
    lines_to_yield = []
    select_by_occupancy = option is None

    atom_lines = ("ATOM", "HETATM")

    # anisou lines are treated specially
    anisou_lines = ("ANISOU",)

    for resnum, atomnames in register.items():
        for atomname, altlocs in atomnames.items():
            if select_by_occupancy:
                # gathers all alternative locations for the atom
                all_lines = []
                for altloc, lines in altlocs.items():
                    all_lines.extend(lines)

                # identifies the highest occupancy combining dictionary
                # and sorting
                new = {}
                for line_number, line in all_lines:
                    if line.startswith(atom_lines):
                        occupancy_number = line[54:60]
                        list_ = new.setdefault(occupancy_number, [])
                        list_.append((line_number, line))

                    # assumes ANISOU succeed the respective ATOM line
                    elif line.startswith(anisou_lines):
                        list_.append((line_number, line))

                # sort keys by occupancy
                keys_ = sorted(
                    new.keys(), key=lambda x: float(x.strip()), reverse=True
                )

                these_atom_lines = new[keys_[0]]

                # always yield the first line
                lines_to_yield.extend(_remove_altloc(these_atom_lines[0:1]))

                del all_lines, new

            # selected by option:
            else:
                if option in altlocs:
                    # selects the option, that's it
                    lines_to_yield.extend(_remove_altloc(altlocs[option]))

                else:
                    # if the option does not exist, add all altlocs
                    for altloc, lines in altlocs.items():
                        lines_to_yield.extend(lines)

    # add comments
    lines_to_yield.extend(others)

    # lines are sorted to the line number so that the output is sorted
    # the same way as in the input PDB
    lines_to_yield.sort(key=lambda x: x[0])

    # the line number is ignored, only the line is yield
    for line_number, line in lines_to_yield:
        yield line


def select_by_occupancy(fhandle, option=None):
    """
    Selects altloc labels for the entire PDB file.

    Parameters
    ----------
    fhandle : an iterable giving the PDB file line-by-line.

    option : str or `None`.
        The alternative location identifier to select. By default
        (`None`) selects the alternative location with highest
        occupancy. In this case, if the different alternative locations
        have the same occupancy, selects the one that comes first.
        Selecting by highest occupancy removes all altloc labels for all
        atoms. Provide an option (e.g. 'A') to select only atoms with
        altloc label `A`. If you select `A` and an atom has conformers
        with altlocs `B` and `C`, both B and C will be kept in the
        output. Despite not an official format, many times alternative
        locations are identified by a blank character ' ' (space), and a
        [A-Z] character.  In these cases, to select the alternative
        location identified by a blank character give `option=' '`.

    Returns
    -------
    generator
        A generator object. To exhaust the generator, that is, to
        process the PDB file (or PDB lines), convert it to a list.

        >>> from pdbtools.pdb_selaltloc import run
        >>> with('input.pdb', 'r') as fin:
        >>>     processed_lines = list(run(fin))

        For more example see:

        >>> import pdbtools
        >>> help(pdbtools)
    """
    records = ("ATOM", "HETATM", "ANISOU")
    terminators = ("TER", "END", "CONECT", "END", "ENDMDL", "MODEL")

    # register atom information
    register = dict()

    # register comment lines
    others = []

    # register current chain
    chain = None
    prev_chain = None

    # keep record of the line number. This will be used to sort lines
    # after selecting the desired alternative location
    nline = 0

    # the loop will collect information on the different atoms
    # throughout the PDB file until a new chain or any terminal line is
    # found. At that point, the collected information is flushed because
    # all altlocs for that block have been defined.
    for line in fhandle:
        nline += 1

        if line.startswith(records):
            # here resnum + insertion code are taken to identify
            # different residues
            resnum = line[22:27]
            atomname = line[12:16]
            altloc = line[16]
            chain = line[21:22]

            # flush lines because we enter a new chain
            if chain != prev_chain:
                # the "yield from" statement is avoided to keep
                # compatibility with Python 2.7
                for _line in _flush(register, option, others):
                    yield _line

                # Python 2.7 compatibility. Do not use .clear() method
                # restart help variables
                del register, others
                register, others = dict(), []

            # organizes information hierarchically
            resnum_d = register.setdefault(resnum, {})
            atomname_d = resnum_d.setdefault(atomname, {})
            altloc_d = atomname_d.setdefault(altloc, [])

            # adds info to dictionary
            altloc_d.append((nline, line))

        # flush information because we reached the end of a block
        elif line.startswith(terminators):
            for _line in _flush(register, option, others):
                yield _line

            del register, others
            register, others = dict(), []

            yield line  # yield the current line after flush

        else:
            # append comments to flush list
            # The reason to add comments to a list instead of yielding
            # them directly is to cover the possibility of having
            # comments in the middle of the PDB file. Obviously is this
            # extremely unlikely. But just in case...
            others.append((nline, line))

        prev_chain = chain

    # at the end of the PDB, flush the remaining lines
    for _line in _flush(register, option, others):
        yield _line


def fetch_updated_cif(pdb_id, cif_fname):
    """
    Fetch updated cif from PDBE database.

    Parameters
    ----------
    pdb_id : str
        PDB ID
    cif_fname : str or Path
        name of the output cif file

    Returns
    -------
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
    return Path(cif_fname)


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


def check_big_uni(ats_dict, uniprot_id):
    """
    Checks if uniprot id has residue IDs > 9999 in the mmcif atom_site dict

    Parameters
    ----------
    ats_dict : dict
        mmcif atom_site dictionary
    uniprot_id : str
        uniprot ID

    Returns
    -------
    big_uni : bool
        True if uniprot ID has residue IDs > 9999
    """
    len_sifts_mapping = len(ats_dict["auth_seq_id"])
    big_uni = False
    for residx in range(len_sifts_mapping):
        curr_uniprot_id = ats_dict["pdbx_sifts_xref_db_acc"][residx]
        if curr_uniprot_id == uniprot_id:
            if int(ats_dict["pdbx_sifts_xref_db_num"][residx]) > 9999:
                big_uni = True
                break
    return big_uni


def convert_cif_to_pdbs(cif_fname, pdb_id, uniprot_id):
    """
    Converts a cif file into a pdb file for each chain matching the uniprot_id

    Parameters
    ----------
    cif_fname : str or Path
        input cif file
    pdb_id : str
        PDB ID
    uniprot_id : str
        uniprot ID to be used

    Returns
    -------
    out_pdb_fnames : list
        list of pdb filenames
    """
    cif_dict = get_cif_dict(cif_fname)
    ats_dict = cif_dict[pdb_id.upper()]["_atom_site"]
    len_sifts_mapping = len(ats_dict["auth_seq_id"])
    # initialising lists
    out_pdb_fnames, out_pdb_lines, atom_ids = [], [], []
    big_uni = check_big_uni(ats_dict, uniprot_id)
    if big_uni:
        log.info(f"uniprot id {uniprot_id} in {pdb_id} has residue IDs > 9999")
    else:
        # iterating over the atomsite_dict
        for residx in range(len_sifts_mapping):
            # extracting key variables from the atom_site dict
            atom_keyword = ats_dict["group_PDB"][residx]
            resid = ats_dict["pdbx_sifts_xref_db_num"][residx]
            curr_uniprot_id = ats_dict["pdbx_sifts_xref_db_acc"][residx]
            atom_symbol = ats_dict["type_symbol"][residx]
            chain = (
                ats_dict["auth_asym_id"][residx]
                if len(ats_dict["auth_asym_id"][residx]) == 1
                else ats_dict["label_asym_id"][residx]
            )
            model_id = int(ats_dict["pdbx_PDB_model_num"][residx])
            # given the valuse of these variables, check if we have to
            # consider this line
            if (
                atom_keyword == "ATOM"
                and resid != "?"
                and curr_uniprot_id == uniprot_id
                and atom_symbol not in ["H", "D"]
                and model_id == 1
            ):
                # getting the correct pdb filename to write on
                pdb_fname = Path(f"{pdb_id}-{chain}.pdb")
                if pdb_fname not in out_pdb_fnames:
                    out_pdb_fnames.append(pdb_fname)
                    atom_ids.append(0)
                # updating atom_id
                atom_idx = out_pdb_fnames.index(pdb_fname)
                atom_ids[atom_idx] += 1
                atom_id = atom_ids[atom_idx]
                # literal data
                atom_name = ats_dict["label_atom_id"][residx]
                alt_id = (
                    ats_dict["label_alt_id"][residx]
                    if ats_dict["label_alt_id"][residx] != "."
                    else " "
                )
                resname = ats_dict["label_comp_id"][residx]
                # ins_code = (
                #     ats_dict["pdbx_PDB_ins_code"][residx]
                #     if ats_dict["pdbx_PDB_ins_code"][residx] != "?"
                #     else " "
                # )
                # numbers
                x = "{:.3f}".format(float(ats_dict["Cartn_x"][residx]))
                y = "{:.3f}".format(float(ats_dict["Cartn_y"][residx]))
                z = "{:.3f}".format(float(ats_dict["Cartn_z"][residx]))
                occ = "{:.2f}".format(float(ats_dict["occupancy"][residx]))
                bfactor = "{:.2f}".format(
                    float(ats_dict["B_iso_or_equiv"][residx])
                )
                # creating new line
                new_line = (
                    f"{atom_keyword}  {atom_id:>5}  {atom_name:<3}{alt_id:<1}"
                    f"{resname:<3} {chain}{resid:>4} "  # removed ins_code
                    f"{x:>11}{y:>8}{z:>8}{occ:>6}{bfactor:>6}"
                    f"{os.linesep}"
                )
                # appending new line to the correct pdb file
                if len(out_pdb_lines) <= atom_idx:
                    out_pdb_lines.append([new_line])
                else:
                    out_pdb_lines[atom_idx].append(new_line)
        # write files
        for pdb_fname, pdb_lines in zip(out_pdb_fnames, out_pdb_lines):
            with open(pdb_fname, "w") as wfile:
                for new_line in pdb_lines:
                    wfile.write(new_line)
    return out_pdb_fnames


def fetch_pdb_files(pdb_to_fetch, uniprot_id):
    """
    Fetches the pdb files from PDBe database.

    Parameters
    ----------
    pdb_to_fetch : list
        list of pdb hits to fetch

    Returns
    -------
    validated_pdbs : list
        list of tuples (pdb_file, cif_file, hit)
    """
    validated_pdb_and_cifs = []
    valid_pdb_set = set()  # set of valid pdb IDs
    for hit in pdb_to_fetch:
        pdb_id = hit["pdb_id"]
        chain_id = hit["chain_id"]
        cif_fname = f"{pdb_id}_updated.cif"
        # if the cif file has not been downloaded yet, download it
        if cif_fname not in os.listdir():
            cif_f = fetch_updated_cif(pdb_id, cif_fname)
            pdb_files = convert_cif_to_pdbs(cif_f, pdb_id, uniprot_id)
            log.info(f"converted cif to pdb files: {pdb_files}")
        else:
            cif_f = Path(cif_fname)
        pdb_fname = f"{pdb_id}-{chain_id}.pdb"
        pdb_f = Path(pdb_fname)
        if pdb_f.exists():
            validated_pdb_and_cifs.append((pdb_f, cif_f, hit))
            if pdb_id not in valid_pdb_set:
                valid_pdb_set.add(pdb_id)
    return validated_pdb_and_cifs


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
    uniprot_id,
    check_pdb,
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
    uniprot_id : str
        Uniprot ID.
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
        chain_id = hit["chain_id"]
        coverage = hit["coverage"]
        resolution = hit["resolution"]
        exp_method = hit["experimental_method"]
        if check_pdb:
            # check coverage value
            if coverage > coverage_cutoff:
                check_list.append(True)
            else:
                check_list.append(False)
                reason = "coverage"
            # check resolution value
            if resolution is None:
                if "NMR" in exp_method:
                    check_list.append(True)
                else:
                    check_list.append(False)
                    reason = "None resolution"
            elif resolution < resolution_cutoff:
                check_list.append(True)
            else:
                check_list.append(False)
                reason = "resolution"

        # check chain ID not longer than 1 character
        # this check holds also if check_pdb is False
        if len(chain_id) == 1:
            check_list.append(True)
        else:
            check_list.append(False)
            reason = "chain ID too big"

        # append pdb to fetch list if all checks passed
        if all(check_list):
            pdbs_to_fetch.append(hit)
        else:
            log.debug(f"{pdb_id}-{chain_id} failed validation ({reason})")
    log.info(f"Found {len(pdbs_to_fetch)} valid PDBs to fetch")
    # downloading a list of good pdbs
    validated_pdbs_and_cifs = fetch_pdb_files(
        pdbs_to_fetch[:max_pdb_num], uniprot_id
    )
    log.info(f"Fetched {len(validated_pdbs_and_cifs)} valid PDBs")
    return validated_pdbs_and_cifs


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
    occ_pdb_f = occ_pdb(pdb_fname)
    tidy_pdb_f = tidy_pdb(occ_pdb_f)

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
    pdb_f : Path
        Path to best PDB file.
    cif_f : Path
        Path to best CIF file.
    hit : dict
        Best hit.
    filtered_interfaces : dict
        Dictionary of filtered interfaces.
    """
    log.info("Selecting pdb retaining the most interfaces")
    cif_f, pdb_f, hit, filtered_interfaces = None, None, None, None
    if validated_pdbs != []:
        max_nint = 0
        for curr_pdb, curr_cif_f, curr_hit in validated_pdbs:
            chain_id = curr_hit["chain_id"]

            # preprocessing pdb file
            tidy_pdb_f = preprocess_pdb(curr_pdb, chain_id)

            try:
                mdu = mda.Universe(tidy_pdb_f)
            except Exception as e:
                log.error(f"Error loading {tidy_pdb_f}: {e}")
                continue
            selection_string = f"name CA and chainID {chain_id.upper()}"
            pdb_resids = mdu.select_atoms(selection_string).resids
            tmp_filtered_interfaces = filter_interfaces(
                interface_residues, pdb_resids
            )
            curr_nint = len(tmp_filtered_interfaces)
            if curr_nint > max_nint:  # update "best" hit
                max_nint = curr_nint
                filtered_interfaces = tmp_filtered_interfaces.copy()
                pdb_f = tidy_pdb_f
                cif_f = curr_cif_f
                hit = curr_hit
        # unlink pdb and cif files
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
    check_pdb = True
    if pdb_to_use:
        pdb_to_use = pdb_to_use.lower()
        check_pdb = False
    if chain_to_use:
        chain_to_use = chain_to_use.upper()
    pdb_list = filter_pdb_list(pdb_dict[uniprot_id], pdb_to_use, chain_to_use)

    validated_pdbs_and_cifs = validate_api_hit(pdb_list, uniprot_id, check_pdb)

    pdb_f, cif_f, top_hit, filtered_interfaces = get_maxint_pdb(
        validated_pdbs_and_cifs,
        interface_residues,
    )

    if pdb_f is None or cif_f is None:
        log.warning(f"Could not fetch PDB/mmcif file for {uniprot_id}")
        return None, None, None

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
