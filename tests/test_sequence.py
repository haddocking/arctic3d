from Bio import SeqIO
from pathlib import Path
import os
import pytest

from arctic3d.modules.sequence import cycle_alignment

from . import golden_data


@pytest.fixture
def inp_fasta():
    return Path(golden_data, "1rypB_r_b.fasta")


def test_cycle_alignment(inp_fasta):
    """Test cycle_alignment function."""
    list_of_seqs = SeqIO.parse(open(inp_fasta), "fasta")
    ref_seq = "MTDRYSFSLTTFSPSGKLGQIDYALTAVKQGVTSLGIKATNGVVIATEKKSSSPLAMSET"
    max_id_chain, max_id = cycle_alignment(list_of_seqs, ref_seq, "aln")
    aln_content = open("aln", "r").read()
    assert max_id_chain == "B"
    assert max_id == 1.0
    first_line = aln_content.split(os.linesep)[0]
    exp_first_line = f"target            0 {ref_seq}"
    assert first_line == exp_first_line
    os.remove("aln")
