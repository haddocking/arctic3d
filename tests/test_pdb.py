from pathlib import Path
import os
import pytest

from arctic3d.modules.pdb import (
    filter_pdb_list,
    get_best_pdb,
    get_maxint_pdb,
    keep_atoms,
    occ_pdb,
    selchain_pdb,
    selmodel_pdb,
    tidy_pdb,
    validate_api_hit,
)

from . import golden_data


@pytest.fixture
def inp_pdb():
    return Path(golden_data, "1rypB_r_b.pdb")


@pytest.fixture
def tricky_pdb():
    return Path(golden_data, "pdb4xoj.pdb")


@pytest.fixture
def inp_pdb_3psg():
    return Path(golden_data, "pdb3psg.pdb")


@pytest.fixture
def inp_cif_3psg():
    return Path(golden_data, "3psg_updated.cif")


@pytest.fixture
def inp_pdb_data():
    return Path(golden_data, "pdb_data_P40202.json")


@pytest.fixture
def pdb_hit_no_resolution():
    hit = {
        "end": 951,
        "chain_id": "A",
        "pdb_id": "2gsx",
        "start": 1,
        "unp_end": 971,
        "coverage": 0.939,
        "unp_start": 21,
        "resolution": None,
        "experimental_method": "X-ray solution scattering",
        "tax_id": 9606,
    }
    return hit


@pytest.fixture
def good_hits():
    hits_list = [
        {
            "end": 246,
            "chain_id": "A",
            "pdb_id": "4xoj",
            "start": 1,
            "unp_end": 246,
            "coverage": 1,
            "unp_start": 1,
            "resolution": 0.91,
            "experimental_method": "X-ray diffraction",
            "tax_id": 9913,
        },
        {
            "end": 246,
            "chain_id": "A",
            "pdb_id": "6sy3",
            "start": 1,
            "unp_end": 246,
            "coverage": 1,
            "unp_start": 1,
            "resolution": 0.95,
            "experimental_method": "X-ray diffraction",
            "tax_id": 9913,
        },
        {
            "end": 246,
            "chain_id": "B",
            "pdb_id": "4gux",
            "start": 1,
            "unp_end": 246,
            "coverage": 1,
            "unp_start": 1,
            "resolution": 1.803,
            "experimental_method": "X-ray diffraction",
            "tax_id": 9913,
        },
        {
            "end": 246,
            "chain_id": "C",
            "pdb_id": "4gux",
            "start": 1,
            "unp_end": 246,
            "coverage": 1,
            "unp_start": 1,
            "resolution": 1.803,
            "experimental_method": "X-ray diffraction",
            "tax_id": 9913,
        },
    ]
    return hits_list


@pytest.fixture
def example_interfaces():
    interfaces = {
        "P01024": [103, 104, 105],
        "P-dummy": [103, 104, 105, 1049, 1050],
    }
    return interfaces


def test_selchain_pdb(inp_pdb):
    pdb = selchain_pdb(inp_pdb, "B")
    assert pdb.exists()
    pdb.unlink()


def test_tidy_pdb(inp_pdb):
    pdb = tidy_pdb(inp_pdb)
    assert pdb.exists()
    pdb.unlink()


def test_occ_pdb(tricky_pdb):
    pdb = occ_pdb(tricky_pdb)
    assert pdb.exists()
    # check that LYS84 is correctly processed
    obs_lys84_lines = []
    with open(pdb, "r") as fh:
        for line in fh:
            if line.startswith("ATOM"):
                if line.split()[5] == "84":
                    obs_lys84_lines.append(line)
    exp_lys84_lines = [
        "ATOM    536  N   LYS A  84      -4.827  -2.055  19.735  1.00  9.68           N  ",  # noqa: E501
        "ATOM    537  CA  LYS A  84      -4.644  -3.431  19.279  0.70 10.05           C  ",  # noqa: E501
        "ATOM    539  C   LYS A  84      -3.278  -3.508  18.640  1.00  9.20           C  ",  # noqa: E501
        "ATOM    540  O   LYS A  84      -2.951  -2.592  17.884  1.00  9.59           O  ",  # noqa: E501
        "ATOM    541  CB  LYS A  84      -5.744  -3.907  18.283  0.70 10.32           C  ",  # noqa: E501
        "ATOM    543  CG  LYS A  84      -7.150  -3.728  18.755  0.70 12.75           C  ",  # noqa: E501
        "ATOM    545  CD  LYS A  84      -8.204  -4.287  17.849  0.70 14.48           C  ",  # noqa: E501
        "ATOM    547  CE  LYS A  84      -9.569  -4.039  18.445  0.70 21.19           C  ",  # noqa: E501
        "ATOM    549  NZ  LYS A  84     -10.614  -4.841  17.764  0.70 22.64           N  ",  # noqa: E501
    ]
    exp_lys84_lines = [ln + os.linesep for ln in exp_lys84_lines]
    assert obs_lys84_lines == exp_lys84_lines
    pdb.unlink()


def test_keep_atoms(inp_pdb):
    pdb = keep_atoms(inp_pdb)
    assert pdb.exists()
    pdb.unlink()


def test_selmodel_pdb(inp_pdb):
    pdb = selmodel_pdb(inp_pdb, "1")
    assert pdb.exists()
    pdb.unlink()


def test_validate_api_hit(pdb_hit_no_resolution):
    """Test validate_api_hit."""
    validated_pdbs = validate_api_hit([pdb_hit_no_resolution], "P20023")
    assert (
        validated_pdbs == []
    )  # this is empty because resolution is None and exp != NMR
    # change resolution to 1.0
    pdb_hit_no_resolution["resolution"] = 1.0
    validated_pdbs = validate_api_hit([pdb_hit_no_resolution], "P20023")
    pdb, cif, dict = validated_pdbs[0]
    assert pdb.name == "2gsx-A.pdb"
    assert cif.name == "2gsx_updated.cif"
    assert dict == pdb_hit_no_resolution


def test_validate_api_hit_nmr(pdb_hit_no_resolution):
    """Test validate_api_hit with NMR data."""
    pdb_hit_no_resolution["experimental_method"] = "Solution NMR"
    # NMR structures have no resolution but should be accepted
    validated_pdbs = validate_api_hit([pdb_hit_no_resolution], "P20023")
    pdb, cif, dict = validated_pdbs[0]
    assert pdb.name == "2gsx-A.pdb"
    assert cif.name == "2gsx_updated.cif"
    assert dict == pdb_hit_no_resolution


def test_get_best_pdb(example_interfaces):
    """Test get_best_pdb."""
    pdb, cif, filtered_interfaces = get_best_pdb("P20023", example_interfaces)
    exp_pdb = Path("P20023-1ghq-B.pdb")
    exp_cif = Path("1ghq_updated.cif")
    exp_interfaces = {"P01024": [103, 104, 105]}
    assert pdb == exp_pdb
    assert cif == exp_cif
    assert filtered_interfaces == exp_interfaces
    exp_pdb.unlink()
    exp_cif.unlink()


def test_get_maxint_pdb_empty():
    """Test get_maxint_pdb with empty output."""
    empty_validated_pdbs = []
    pdb_f, cif_f, top_hit, filtered_interfaces = get_maxint_pdb(
        empty_validated_pdbs, {}
    )
    assert pdb_f is None
    assert cif_f is None
    assert top_hit is None
    assert filtered_interfaces is None


def test_get_maxint_pdb(good_hits, example_interfaces):
    """Test get_maxint_pdb with implicit pdb numbering."""
    validated_pdbs = validate_api_hit(good_hits, "P00760")
    pdb_f, cif_f, top_hit, filtered_interfaces = get_maxint_pdb(
        validated_pdbs, example_interfaces
    )
    assert pdb_f.name == "4xoj-A-occ-tidy.pdb"
    assert cif_f.name == "4xoj_updated.cif"
    assert top_hit["pdb_id"] == "4xoj"
    assert top_hit["chain_id"] == "A"
    assert filtered_interfaces == {"P01024": [103, 104, 105]}


def test_filter_pdb_list(good_hits):
    """Test filter_pdb_list."""
    observed_red_list = filter_pdb_list(good_hits, pdb_to_use="1abc")
    expected_red_list = []
    assert observed_red_list == expected_red_list
    observed_red_list = filter_pdb_list(good_hits, pdb_to_use="6sy3")
    expected_red_list = [good_hits[1]]
    assert observed_red_list == expected_red_list
    # testing chain information
    observed_red_list = filter_pdb_list(
        good_hits, pdb_to_use="4gux", chain_to_use="C"
    )
    expected_red_list = [good_hits[3]]
    assert observed_red_list == expected_red_list


def test_pdb_data(inp_pdb_data):
    """Test pdb_data input json file."""
    orig_interfaces = {"P00441": [85, 137, 138]}
    pdb, cif, filtered_interfaces = get_best_pdb(
        "P40202", orig_interfaces, pdb_data=inp_pdb_data
    )

    assert filtered_interfaces == orig_interfaces
    pdb.unlink()
    cif.unlink()
