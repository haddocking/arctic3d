import tempfile

from Bio import Align
from Bio.Align import substitution_matrices
from pdbtools.pdb_tofasta import pdb_to_fasta


def load_seq(fasta_file):
    """
    Load sequence.

    Parameters
    ----------
    fasta_file : str
        Fasta file.

    Returns
    -------
    fasta_seq : str
        Fasta sequence.

    """
    pass


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
        fasta_fh = tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False)
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
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM80")
    alns = aligner.align(seq1, seq2)
    top_aln = alns[0]
    with open(aln_fname, "w") as fh:
        fh.write(str(top_aln))
    return aln_fname, top_aln
