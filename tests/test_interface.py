from pathlib import Path

import pytest

from arctic3d.modules.interface import parse_out_uniprot, read_interface_residues

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


def test_parse_out_uniprot():
    uniprot_strings = [None, "P00760", "P00760,P00974"]
    expected_uniprot_strings = [[], ["P00760"], ["P00760", "P00974"]]
    observed_uniprot_strings = []
    for string in uniprot_strings:
        obs_list = parse_out_uniprot(string)
        observed_uniprot_strings.append(obs_list)
    assert expected_uniprot_strings == expected_uniprot_strings
