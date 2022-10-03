import os
from pathlib import Path

import pytest

from arctic3d.modules.pdb import (
    fetch_pdbrenum,
    filter_pdb_list,
    get_best_pdb,
    get_maxint_pdb,
    keep_atoms,
    occ_pdb,
    output_pdb,
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
    ]
    return hits_list


def test_fetch_pdbrenum():
    pdb = fetch_pdbrenum("1crn")
    assert pdb.exists()
    pdb.unlink()


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
    pdb, filtered_interfaces = get_best_pdb("P20023", {"P01024": [103, 104, 105]})
    assert pdb is None
    assert filtered_interfaces is None


def test_get_maxint_pdb():
    """Test get_maxint_pdb."""
    empty_validated_pdbs = []
    pdb_f, top_hit, filtered_interfaces = get_maxint_pdb(empty_validated_pdbs, {})
    assert pdb_f is None
    assert top_hit is None
    assert filtered_interfaces is None

    # TODO: test the non-empty case as well


def test_filter_pdb_list(good_hits):
    """Test filter_pdb_list."""
    observed_red_list = filter_pdb_list(good_hits, "1abc")
    expected_red_list = []
    assert observed_red_list == expected_red_list
    observed_red_list = filter_pdb_list(good_hits, "6sy3")
    expected_red_list = [good_hits[1]]
    assert observed_red_list == expected_red_list


def test_output_pdb(inp_pdb):
    """Test output_pdb."""
    example_res_probs = {1: {3: 0.2, 4: 0.75}}
    output_files = output_pdb(inp_pdb, example_res_probs)
    original_content = open(inp_pdb, "r").read().split(os.linesep)
    original_content = [el for el in original_content if el.startswith("ATOM")]
    # check file existence
    assert Path.exists(output_files[0])
    observed_content = open(output_files[0], "r").read().split(os.linesep)
    observed_content = [el for el in observed_content if el.startswith("ATOM")]
    # assert equal length
    assert len(original_content) == len(observed_content)
    # assert equal content (except for b factors)
    for ln_id in range(len(observed_content)):
        assert original_content[ln_id][:60] == observed_content[ln_id][:60]
    os.unlink(output_files[0])
