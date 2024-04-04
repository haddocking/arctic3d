import logging
import tempfile

from Bio import Align
from Bio.Align import substitution_matrices
from pdbtools.pdb_tofasta import pdb_to_fasta

import os
import shutil

log = logging.getLogger("arctic3d.log")

LETTERS = [
    "A",
    "B",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
    "Z",
]


def to_fasta(pdb_f, temp):
    """
    Convert PDB to Fasta.

    Parameters
    ----------
    pdb_f : str
        PDB filename.

    Returns
    -------
    fasta_fh : str
        Fasta sequence.

    """
    if temp:
        fasta_fh = tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        )
    else:
        fasta_fh = open(pdb_f.with_suffix(".fasta"), "w")

    with open(pdb_f, "r") as inp_fh:
        for line in pdb_to_fasta(inp_fh, multi=False):
            fasta_fh.write(line)

    fasta_fh.close()

    return fasta_fh


def align_sequences(seq1, seq2):
    """
    Performs a pairwise alignment between two sequences.

    Parameters
    ----------
    seq1 : str
        sequence
    seq2 : str
        sequence

    Returns
    -------
    aln_fname : str
        alignment file name
    top_aln : str
        top alignment
    """
    aln_fname = "blosum80.aln"
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -1.01
    aligner.extend_gap_score = -1.000
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

    alns = aligner.align(seq1, seq2)
    top_aln = alns[0]
    with open(aln_fname, "w") as fh:
        fh.write(str(top_aln))
    return aln_fname, top_aln


def cycle_alignment(fasta_sequences, ref_seq, output_aln_fname):
    """
    Perform pairwise alignment between a reference sequence and a list of
     sequences.

    Parameters
    ----------
    fasta_sequences : list
        List of sequences

    ref_seq : str
        Reference sequence

    output_aln_fname : str
        Output alignment file name
    """
    # looping over sequences
    max_id = -1.0
    for fasta in fasta_sequences:
        name, seq = fasta.id, str(fasta.seq)
        try:
            aln_fname, top_aln = align_sequences(ref_seq, seq)
        except Exception as e:
            log.warning(
                f"Error aligning sequence {name} to reference."
                "Is it DNA/RNA? Skipping alingment..."
                f"Error: {e}"
            )
            identity = -1.0
            continue
        identity = str(top_aln).count("|") / float(min(len(ref_seq), len(seq)))
        # compute percentage and logging
        perc_identity = identity * 100
        log.info(f"sequence {name} has {perc_identity:.2f}% sequence identity")
        if identity > max_id:
            max_id = identity
            max_id_chain = name.split("|")[1]
            shutil.copy(aln_fname, output_aln_fname)
    os.unlink(aln_fname)
    return max_id_chain, max_id


def extract_aln_string(pdb_numb_ln, nlines):
    """
    Extracts sequence strings from the alignment file.

    Parameters
    ----------
    pdb_numb_ln : list
        list of alignment lines

    nlines : list
        list of numbers of lines tp be extracted

    Returns
    -------
    aln_lines : list
        list of alignment lines
    """
    aln_lines = []
    for n_ln in nlines:
        splt_ln = pdb_numb_ln[n_ln].split()
        if len(splt_ln) > 2:
            aln_lines.append(splt_ln[2])
    return aln_lines
