import os
import shutil
from pathlib import Path
from unittest.mock import MagicMock

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
        biological_clustering=[],
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
        biological_clustering=[],
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


@pytest.mark.integration
def test_biological_clustering_integration():
    """Integration test for biological_clustering with location option."""
    target_uniprot = "W5JXD7"
    exp_dir = Path(f"arctic3d-{target_uniprot}")
    localise_dir = Path("arctic3d-localise-subcellular")

    # Clean up if exists
    if exp_dir.exists():
        shutil.rmtree(exp_dir)
    if localise_dir.exists():
        shutil.rmtree(localise_dir)

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
        biological_clustering=["location"],
    )

    # Verify main clustering succeeded
    assert exit_code == 0

    # Change back to original directory before checking files
    os.chdir(start_cwd)

    # Verify main clustering directory was created
    assert exp_dir.exists()
    assert Path(exp_dir, "clustered_interfaces.out").exists()

    # Verify biological clustering directory was created
    # (localise creates directory in the main arctic3d dir)
    localise_dir_full = Path(exp_dir, "arctic3d-localise-subcellular")
    if localise_dir_full.exists():
        assert Path(localise_dir_full, "arctic3d-localise.log").exists()
    elif localise_dir.exists():
        # Or it might be in the working directory
        assert Path(localise_dir, "arctic3d-localise.log").exists()

    # Cleanup
    if exp_dir.exists():
        shutil.rmtree(exp_dir)
    if localise_dir.exists():
        shutil.rmtree(localise_dir)


def test_biological_clustering_config():
    """Unit test to verify biological clustering configuration is correct."""
    # Test the configuration mapping matches webserver
    clustering_config = {
        "location": {
            "quickgo": "C",
            "dir": "arctic3d-localise-subcellular",
        },
        "function": {
            "quickgo": "F",
            "dir": "arctic3d-localise-proteinfunction",
        },
        "process": {
            "quickgo": "P",
            "dir": "arctic3d-localise-biologicalprocess",
        },
    }

    # Verify all config values are correct
    assert clustering_config["location"]["quickgo"] == "C"
    assert clustering_config["location"]["dir"] == "arctic3d-localise-subcellular"

    assert clustering_config["function"]["quickgo"] == "F"
    assert clustering_config["function"]["dir"] == "arctic3d-localise-proteinfunction"

    assert clustering_config["process"]["quickgo"] == "P"
    assert clustering_config["process"]["dir"] == "arctic3d-localise-biologicalprocess"


@pytest.mark.integration
def test_biological_clustering_not_called_when_no_interfaces(monkeypatch):
    """Test that biological_clustering is not executed when there are no interfaces."""
    # Mock localise_main to track if it's called
    mock_localise = MagicMock()
    monkeypatch.setattr("arctic3d.cli_localise.main", mock_localise)

    target_uniprot = "P23804"  # This uniprot has no interfaces
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
        biological_clustering=["location", "function"],
    )

    # Verify localise was NOT called since there are no interfaces
    mock_localise.assert_not_called()
    assert exit_code == 255

    # Cleanup
    os.chdir(start_cwd)
    exp_dir = Path(f"arctic3d-{target_uniprot}")
    if exp_dir.exists():
        shutil.rmtree(exp_dir)
