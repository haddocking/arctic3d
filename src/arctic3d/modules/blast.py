"""Function to BLAST input sequence and return accession id."""

import logging
import os
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path

from Bio.Blast import NCBIWWW
from lxml import etree as ET

log = logging.getLogger("arctic3d.log")


def run_blast(fasta_f, db=None):
    """
    Run BLAST.

    Parameters
    ----------
    local : bool
        Run BLAST locally.

    """
    if db:
        log.info(f"Running BLAST locally against {Path(db).name}")
        accession_id = blast_local(fasta_f, db)
    else:
        log.info("Running BLAST remotely...")
        accession_id = blast_remote(fasta_f)

    return accession_id.split(".")[0]


def get_blast_exec():
    """
    Get BLAST executable.

    Returns
    -------
    blastp_exec : str
        BLASTp executable.

    """
    if shutil.which("blastp"):
        blastp_exec = "blastp"
    else:
        ncbi_blast_path = [
            f
            for f in Path(__file__).parent.parent.parent.glob("ncbi-blast*")
            if f.is_dir()
        ][0]
        blastp_exec = Path(ncbi_blast_path, "bin/blastp")
        if not blastp_exec.exists():
            log.error("Could not find blastp executable")
    return blastp_exec


def blast_local(fasta_file, db):
    """
    Blast sequence against Uniprot locally.

    Parameters
    ----------
    fasta_file : str
        Fasta filename

    Returns
    -------
    uniprot_id : str
        Uniprot ID.

    """
    blastp_exec = get_blast_exec()
    cmd = f"{blastp_exec} -query {fasta_file} -db {db} -outfmt 6"

    p = subprocess.run(
        shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    if p.returncode != 0:
        log.error(p.stderr.decode())
        raise Exception("BLAST failed")

    out = p.stdout.decode("utf-8").split(os.linesep)
    uniprot_id = out[0].split("\t")[1]
    return uniprot_id


def blast_remote(fasta_file: str) -> str:
    """
    Blast sequence.

    Parameters
    ----------
    fasta_seq : str
        Fasta filename.

    Returns
    -------
    uniprot_id : str
        Uniprot ID.

    """
    blast_res_handle = NCBIWWW.qblast(
        "blastp", "swissprot", fasta_file, hitlist_size=50
    )

    # TODO: Handle scenario in which the `qblast` call fails

    with tempfile.NamedTemporaryFile(
        mode="w+", delete=True, suffix=".xml"
    ) as temp:
        blast_res = blast_res_handle.read()
        temp.write(blast_res)
        temp.flush()
        accession_id = parse_xml(temp.name)

    return accession_id


def parse_xml(xml_file: str) -> str:
    """Parse the BLAST XML file and return the first (?) accession ID."""
    tree = ET.parse(source=xml_file, parser=ET.XMLParser(encoding="utf-8"))
    root = tree.getroot()

    # root [BlastOutput_iterations] [Iteration] [Iteration_hits] \
    #   [Hit #2] [Hit_accession]
    # using second hit as the first is the input
    # instead of Hit_accession, [1] for [Hit_id] can be used
    accession_id = root[8][0][4][1][3].text

    return accession_id
