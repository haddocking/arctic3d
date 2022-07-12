"""Function to BLAST input sequence and return accession id."""
import os

from Bio.Blast import NCBIWWW
from defusedxml import lxml as ET


def blast_seq(fasta_seq):
    """
    Blast sequence.

    Parameters
    ----------
    fasta_seq : str
        Fasta sequence.

    Returns
    -------
    uniprot_id : str
        Uniprot ID.

    """
    blast_res_handle = NCBIWWW.qblast("blastp", "nr", fasta_seq, hitlist_size=50)

    # temp file for storing results
    with open("blast_res.xml", "w") as save_output:
        blast_res = blast_res_handle.read()
        save_output.write(blast_res)

    tree = ET.parse("blast_res.xml")
    root = tree.getroot()

    # root [BlastOutput_iterations] [Iteration] [Iteration_hits] [Hit #2] [Hit_accession]
    # using second hit as the first is the input
    # instead of Hit_accession, [1] for [Hit_id] can be used
    accession_id = root[8][0][4][1][3]

    # cleaning of temporary result file
    os.remove("blast_res.xml")

    return accession_id
