from pathlib import Path

import pytest

from arctic3d.modules.interface import read_interface_residues

from . import golden_data


def test_read_int_file_nonexisting():
    """Test error on non-existing path."""
    non_ex_path = "../dummy"
    with pytest.raises(Exception):
        read_interface_residues(non_ex_path)


def test_read_int_file():
    int_file_path = Path(golden_data, "interface_file.txt")
    obs_interface_dict = read_interface_residues(int_file_path)
    exp_interface_dict = {
        "int_1": [16, 18, 19],
        "int_2": [21, 22, 24, 26, 27, 29, 30, 31, 36, 38, 39, 40],
        "int_3": [16, 17, 18],
        "int_4": [20, 24],
        "int_5": [1, 2, 3],
        "int_6": [1, 2],
    }
    assert obs_interface_dict == exp_interface_dict
