import tempfile

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
