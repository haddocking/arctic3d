import os
import shutil
from pathlib import Path


def setup_output_folder(uniprot_id, input_files):
    """Sets up output folder"""
    if uniprot_id is None:
        key = "custom"
    else:
        key = uniprot_id
    run_dir = Path(f"arctic3d-{key}")
    if os.path.exists(run_dir):
        raise Exception(f"{run_dir} already exists!")
    os.mkdir(run_dir)
    datadir = Path(run_dir, "input_data")
    os.mkdir(datadir)
    for filepath in input_files:
        filename = filepath.name
        shutil.copy(filepath, Path(datadir, filename))
    os.chdir(run_dir)
