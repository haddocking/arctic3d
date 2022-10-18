from pathlib import Path

import pytest

from arctic3d.modules.pdb import (
    filter_pdb_list,
    get_best_pdb,
    get_maxint_pdb,
    keep_atoms,
    occ_pdb,
    selchain_pdb,
    tidy_pdb,
    validate_api_hit,
)

from . import golden_data


@pytest.fixture
def inp_pdb():
    return Path(golden_data, "1rypB_r_b.pdb")


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


def test_selchain_pdb(inp_pdb):
    pdb = selchain_pdb(inp_pdb, "B")
    assert pdb.exists()
    pdb.unlink()


def test_tidy_pdb(inp_pdb):
    pdb = tidy_pdb(inp_pdb)
    assert pdb.exists()
    pdb.unlink()


def test_occ_pdb(inp_pdb):
    pdb = occ_pdb(inp_pdb)
    assert pdb.exists()
    pdb.unlink()


def test_keep_atoms(inp_pdb):
    pdb = keep_atoms(inp_pdb)
    assert pdb.exists()
    pdb.unlink()


def test_validate_api_hit(pdb_hit_no_resolution):
    validated_pdbs = validate_api_hit([pdb_hit_no_resolution])
    assert validated_pdbs == []

    pdb_hit_no_resolution["resolution"] = 1.0
    validated_pdbs = validate_api_hit([pdb_hit_no_resolution])
    pdb, dict = validated_pdbs[0]
    assert pdb.name == "2gsx.pdb"
    assert dict == pdb_hit_no_resolution


def test_get_best_pdb():
    orig_interfaces = {
        "P01024": [103, 104, 105],
        "P-dummy": [103, 104, 105, 1049, 1050],
    }
    pdb, filtered_interfaces = get_best_pdb("P20023", orig_interfaces)
    exp_pdb = Path("P20023-1ghq-B.pdb")
    exp_interfaces = {"P01024": [103, 104, 105]}
    assert pdb == exp_pdb
    assert filtered_interfaces == exp_interfaces
    exp_pdb.unlink()


def test_get_maxint_pdb():
    """Test get_maxint_pdb."""
    empty_validated_pdbs = []
    pdb_f, top_hit, filtered_interfaces = get_maxint_pdb(
        empty_validated_pdbs, {}, uniprot_id=None
    )
    assert pdb_f is None
    assert top_hit is None
    assert filtered_interfaces is None

    # TODO: test the non-empty case as well


def test_filter_pdb_list(good_hits):
    """Test filter_pdb_list."""
    observed_red_list = filter_pdb_list(good_hits, pdb_to_use="1abc")
    expected_red_list = []
    assert observed_red_list == expected_red_list
    observed_red_list = filter_pdb_list(good_hits, pdb_to_use="6sy3")
    expected_red_list = [good_hits[1]]
    assert observed_red_list == expected_red_list
    # testing chain information
    observed_red_list = filter_pdb_list(good_hits, pdb_to_use="4gux", chain_to_use="C")
    expected_red_list = [good_hits[3]]
    assert observed_red_list == expected_red_list
