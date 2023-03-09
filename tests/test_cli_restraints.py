from arctic3d.cli_restraints import (
    main,
    filter_residues_probs,
    generate_restraints,
    compress_tbl_files,
)

from pathlib import Path
import os
import shutil
import glob
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


def test_main():
    """Test main."""
    start_cwd = os.getcwd()
    r1, r2 = Path(golden_data), Path(golden_data)
    run_dir = "arctic3d-restraints"
    main(r1, r2, None, None, run_dir=run_dir, prob_threshold=0.7)
    # check if the zipped tbl files exist
    assert Path("ambig.tbl.tgz").exists()
    # check the correct number of tbl files exist
    ls_tbl = len(glob.glob("ambig*tbl"))
    assert ls_tbl == 4
    os.chdir(start_cwd)
    shutil.rmtree(run_dir)
