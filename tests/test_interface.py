from pathlib import Path

import pytest

from arctic3d.modules.interface import (
    get_interface_residues,
    parse_interface_line,
    parse_out_pdb,
    parse_out_uniprot,
    read_interface_residues,
)

from . import golden_data


@pytest.fixture
def inp_interface_data():
    return Path(golden_data, "interface_data_P40202.json")


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


def test_parse_interface_line():
    """Test parse_interface_line function."""
    interface_lines = ["P00767 1,2,3", "P00767", "P00767 1-3,4"]
    # first string is correct
    exp_interface = "P00767", [1, 2, 3]
    obs_interface = parse_interface_line(interface_lines[0], 0)
    assert exp_interface == obs_interface
    # the other two should throw an exception
    with pytest.raises(Exception):
        parse_interface_line(interface_lines[1], 1)
    with pytest.raises(Exception):
        parse_interface_line(interface_lines[2], 2)


def test_parse_out_uniprot():
    uniprot_strings = [None, "P00760", "P00760,P00974"]
    expected_uniprot_strings = [set([]), set(["P00760"]), set(["P00760", "P00974"])]
    observed_uniprot_strings = []
    for string in uniprot_strings:
        obs_list = parse_out_uniprot(string)
        observed_uniprot_strings.append(obs_list)
    assert expected_uniprot_strings == observed_uniprot_strings


def test_error_parse_out_uniprot():
    out_uniprot_strings = ["P00760+", "P00760/"]
    for string in out_uniprot_strings:
        with pytest.raises(Exception):
            parse_out_uniprot(string)


def test_parse_out_pdb():
    """Test parse_out_pdb function."""
    out_pdb_strings = [None, "1abc", "4B2C,4B1T,4B2A", "4B2C,4B1T,4B2A,4b2a"]
    expected_pdb_sets = [
        set([]),
        set(["1abc"]),
        set(["4b2a", "4b2c", "4b1t"]),
        set(["4b2a", "4b2c", "4b1t"]),
    ]
    observed_pdb_sets = []
    for string in out_pdb_strings:
        obs_pdb_set = parse_out_pdb(string)
        observed_pdb_sets.append(obs_pdb_set)
    assert expected_pdb_sets == observed_pdb_sets


def test_error_parse_out_pdb():
    """Test correct exception of parse_out_pdb function."""
    out_pdb_strings = ["1ab", "1abcd", "4B2C,4B1T,4B2+"]
    for string in out_pdb_strings:
        with pytest.raises(Exception):
            parse_out_pdb(string)


def test_interface_data(inp_interface_data):
    """Test interface_data input json file."""
    obs_interface_residues = get_interface_residues(
        "P40202", None, None, full=False, interface_data=inp_interface_data, ligand="no"
    )
    exp_interface_residues = {
        "P00441": [85, 137, 138, 187, 217, 218, 222, 229, 231, 232],
        "P00445": [136, 137, 138, 187, 217, 218, 226, 229, 230],
        "P40202": [136, 137, 138, 183, 184, 186, 187, 217, 218],
    }
    assert obs_interface_residues == exp_interface_residues


def test_full_interface(inp_interface_data):
    """Test interface_data input json file with full option."""
    obs_interface_residues = get_interface_residues(
        "P40202", None, None, full=True, interface_data=inp_interface_data, ligand="no"
    )
    exp_interface_residues = {
        "P00441-5u9m-A": [85, 137, 138, 187, 217, 218, 222, 229, 231, 232],
        "P00441-5u9m-C": [85, 137, 138, 187, 217, 218, 222, 229, 231, 232],
        "P00445-1jk9-A": [136, 137, 138, 187, 217, 218, 226, 229, 230],
        "P00445-1jk9-C": [136, 137, 138, 187, 217, 218, 226, 229, 230],
        "P40202-1qup-A": [136, 137, 138, 183, 184, 186, 187, 217, 218],
        "P40202-1qup-B": [136, 137, 138, 183, 184, 186, 187, 217, 218],
    }
    assert obs_interface_residues == exp_interface_residues


def test_ligandyes_interface():
    """Test get_interface_residues when ligand == yes."""
    obs_interface_residues = get_interface_residues(
        "P40202", None, None, full=True, interface_data=None, ligand="yes"
    )
    exp_interface_residues = {
        "SO4-1qup-A": [76, 77, 188, 217],
        "SO4-1jk9-D": [226, 229, 230],
        "ZN-1jk9-B": [16],
        "ZN-5u9m-B": [17, 20],
        "CA-1ej8-A": [124, 168],
    }
    assert obs_interface_residues == exp_interface_residues


def test_ligandboth_interface():
    """Test get_interface_residues when ligand == both."""
    obs_interface_residues = get_interface_residues(
        "P40202", None, None, full=True, interface_data=None, ligand="both"
    )
    exp_interface_residues = {
        "P00441-5u9m-A": [85, 137, 138, 187, 217, 218, 222, 229, 231, 232],
        "P00441-5u9m-C": [85, 137, 138, 187, 217, 218, 222, 229, 231, 232],
        "P00445-1jk9-A": [136, 137, 138, 187, 217, 218, 226, 229, 230],
        "P00445-1jk9-C": [136, 137, 138, 187, 217, 218, 226, 229, 230],
        "P40202-1qup-A": [136, 137, 138, 183, 184, 186, 187, 217, 218],
        "P40202-1qup-B": [136, 137, 138, 183, 184, 186, 187, 217, 218],
        "SO4-1qup-A": [76, 77, 188, 217],
        "SO4-1jk9-D": [226, 229, 230],
        "ZN-1jk9-B": [16],
        "ZN-5u9m-B": [17, 20],
        "CA-1ej8-A": [124, 168],
    }
    assert obs_interface_residues == exp_interface_residues
