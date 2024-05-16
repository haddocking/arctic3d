import os
import shutil
from pathlib import Path

import pytest

from arctic3d.cli import main


@pytest.mark.integration
def test_cli_empty():
    """Test main cli with uniprot ID with no interfaces."""
    target_uniprot = "P23804"
    start_cwd = os.getcwd()
    exit_code = main(
        input_arg=target_uniprot,
        db=None,
        interface_file=None,
        out_partner=None,
        out_pdb=None,
        pdb_to_use=None,
        chain_to_use=None,
        run_dir=None,
        interface_data=None,
        pdb_data=None,
        full=None,
        ligand=None,
        linkage_strategy=None,
        threshold=None,
        int_cov_cutoff=None,
        min_clust_size=None,
    )
    # assert exit code
    assert exit_code == 255
    os.chdir(start_cwd)
    exp_dir = Path(f"arctic3d-{target_uniprot}")
    assert exp_dir.exists() is True
    # Check that the log file has been created
    assert Path(exp_dir, "arctic3d.log").exists()
    # remove folder
    if exp_dir.exists():
        shutil.rmtree(exp_dir)


@pytest.mark.integration
def test_cli_full():
    """Test main cli with uniprot ID with one interface."""
    target_uniprot = "W5JXD7"
    exp_dir = Path(f"arctic3d-{target_uniprot}")
    # delete folder if exists
    if exp_dir.exists():
        shutil.rmtree(exp_dir)
    start_cwd = os.getcwd()
    exit_code = main(
        input_arg=target_uniprot,
        db=None,
        interface_file=None,
        out_partner=None,
        out_pdb=None,
        pdb_to_use="3wqb",
        chain_to_use=None,
        run_dir=None,
        interface_data=None,
        pdb_data=None,
        full=None,
        ligand="no",
        linkage_strategy=None,
        threshold=None,
        min_clust_size=1,
        int_cov_cutoff=0.7,
    )
    assert exit_code == 0
    os.chdir(start_cwd)
    assert exp_dir.exists() is True
    # Check that the log file has been created
    assert Path(exp_dir, "arctic3d.log").exists()
    # check content of the clustered interfaces file
    assert Path(exp_dir, "clustered_interfaces.out").exists()
    # check content of the clustered interfaces file
    obs_content = Path(exp_dir, "clustered_interfaces.out").read_text()
    exp_content = f"Cluster 1 -> Q9L5A4-3wqb-A{os.linesep}"
    assert exp_content == obs_content
    # remove folder
    if exp_dir.exists():
        shutil.rmtree(exp_dir)
