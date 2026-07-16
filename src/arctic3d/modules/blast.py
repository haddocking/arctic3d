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
    # qblast expects the query sequence itself, not a file path
    with open(fasta_file) as fh:
        fasta_seq = fh.read()

    blast_res_handle = NCBIWWW.qblast(
        "blastp", "swissprot", fasta_seq, hitlist_size=50
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
    """Parse the BLAST XML file and return the best-hit accession ID."""
    tree = ET.parse(source=xml_file, parser=ET.XMLParser(encoding="utf-8"))
    root = tree.getroot()

    # navigate by tag name (robust to changes in the number of hits and to
    # the exact BlastOutput layout) and return the top hit's accession
    hits = root.findall(".//Hit")
    if not hits:
        raise ValueError(f"No BLAST hits found in {xml_file}")

    accession_id = hits[0].findtext("Hit_accession")
    if not accession_id:
        raise ValueError(f"Could not parse Hit_accession from {xml_file}")

    return accession_id
