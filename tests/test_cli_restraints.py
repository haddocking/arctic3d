import glob
import os
import shutil
from pathlib import Path

import pytest

from arctic3d.cli_restraints import (
    compress_tbl_files,
    filter_residues_probs,
    generate_restraints,
    main,
)

from . import golden_data


def test_filter_residues_probs():
    """Test filter_residues_probs."""
    example_residues_probs = {1: {1: 1.0, 2: 0.8, 3: 0.6}}
    obs_filtered_residues = filter_residues_probs(example_residues_probs, 0.7)
    exp_filtered_residues = {1: [1, 2]}
    assert obs_filtered_residues == exp_filtered_residues


def test_generate_restraints():
    """Test generate_restraints."""
    residues1 = [1, 2]
    residues2 = [5]
    ch1 = "A"
    ch2 = "B"
    ambig_fname = "ambig.tbl"
    generate_restraints(residues1, residues2, ch1, ch2, ambig_fname)
    with open(ambig_fname, "r") as f:
        obs_ambig = f.read()
    exp_ambig = [
        "!HADDOCK AIR restraints for 1st partner",
        "!",
        "assign ( resid 1 and segid A)",
        "       (",
        "        ( resid 5 and segid B)",
        "       )  2.0 2.0 0.0",
        "!",
        "assign ( resid 2 and segid A)",
        "       (",
        "        ( resid 5 and segid B)",
        "       )  2.0 2.0 0.0",
        "!",
        "!",
        "!HADDOCK AIR restraints for 2nd partner",
        "!",
        "assign ( resid 5 and segid B)",
        "       (",
        "        ( resid 1 and segid A)",
        "     or",
        "        ( resid 2 and segid A)",
        "       )  2.0 2.0 0.0",
        "!",
    ]
    exp_ambig = os.linesep.join(exp_ambig)
    exp_ambig += os.linesep
    assert obs_ambig == exp_ambig
    os.unlink(ambig_fname)


def test_compress_tbl_files():
    """Test compress_tbl_files."""
    ambig_fname = "ambig.tbl"
    generate_restraints([1], [2], "A", "B", ambig_fname)
    out_tgz = "ambig.tbl.tgz"
    compress_tbl_files([ambig_fname], out_tgz=out_tgz)
    assert Path(out_tgz).exists()
    os.unlink(ambig_fname)
    os.unlink(out_tgz)


@pytest.mark.integration
def test_main():
    """Test main."""
    start_cwd = os.getcwd()
    r1, r2 = Path(golden_data), Path(golden_data)
    run_dir = "arctic3d-restraints"
    main(r1, r2, None, None, run_dir=run_dir, prob_threshold=0.7)
    # check if the zipped tbl files exist
    assert Path("ambig.tbl.tgz").exists()
    # check if log file exists
    assert Path("arctic3d-restraints.log").exists()
    # check the correct number of tbl files exist
    ls_tbl = len(glob.glob("ambig*tbl"))
    assert ls_tbl == 4
    os.chdir(start_cwd)
    shutil.rmtree(run_dir)
