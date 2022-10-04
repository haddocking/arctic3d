import pytest
from pathlib import Path
import os
from . import golden_data

from arctic3d.modules.output import output_pdb


@pytest.fixture
def inp_pdb():
    return Path(golden_data, "1rypB_r_b.pdb")


def test_output_pdb(inp_pdb):
    """Test output_pdb."""
    example_res_probs = {1: {3: 0.2, 4: 0.75}}
    output_files = output_pdb(inp_pdb, example_res_probs)
    original_content = open(inp_pdb, "r").read().split(os.linesep)
    original_content = [el for el in original_content if el.startswith("ATOM")]
    # check file existence
    assert Path.exists(output_files[0])
    observed_content = open(output_files[0], "r").read().split(os.linesep)
    observed_content = [el for el in observed_content if el.startswith("ATOM")]
    # assert equal length
    assert len(original_content) == len(observed_content)
    # assert equal content (except for b factors)
    for ln_id in range(len(observed_content)):
        assert original_content[ln_id][:60] == observed_content[ln_id][:60]
    os.unlink(output_files[0])
