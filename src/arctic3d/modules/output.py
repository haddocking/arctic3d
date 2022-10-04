import os
import shutil
from pathlib import Path


def setup_output_folder(uniprot_id, input_files):
    """Sets up output folder

    Parameters
    ----------
    uniprot_id : string or None
        uniprot_id of the run
    input_files : dict of Paths
        dict of input files
    Returns
    -------
    copied_input_files : dict of Paths
        dict of copied input files
    """
    if uniprot_id is None:
        key = "custom"
    else:
        key = uniprot_id
    run_dir = Path(f"arctic3d-{key}")
    if os.path.exists(run_dir):
        raise Exception(f"{run_dir} already exists!")
    # setting up the directory
    os.mkdir(run_dir)
    datadir = Path(run_dir, "input_data")
    os.mkdir(datadir)
    for key in input_files:
        filename = input_files[key].name
        filepath = Path(datadir, filename)
        shutil.copy(input_files[key], filepath)
    os.chdir(run_dir)
    # retrieving info about copied input files
    copied_input_files = {}
    for key in input_files:
        filename = input_files[key].name
        copied_input_files[key] = Path("input_data", filename)
    return copied_input_files
