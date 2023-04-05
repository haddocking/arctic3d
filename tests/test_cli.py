from arctic3d.cli import main
import shutil
from pathlib import Path
import os


def test_cli_empty():
    """Test main cli with uniprot ID with no interfaces."""
    target_uniprot = "P23804"
    start_cwd = os.getcwd()
    main(
        target_uniprot,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    )
    os.chdir(start_cwd)
    exp_dir = Path(f"arctic3d-{target_uniprot}")
    assert exp_dir.exists() is True
    # Check that the log file has been created
    assert Path(exp_dir, "arctic3d.log").exists()
    # remove folder
    if exp_dir.exists():
        shutil.rmtree(exp_dir)
