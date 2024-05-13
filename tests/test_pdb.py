import copy
import json
import os
from pathlib import Path

import pytest

from arctic3d.modules.pdb import (
    convert_cif_to_pdbs,
    fetch_pdb_files,
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
def fetch_pdb_files_output():
    return [
        (
            Path(golden_data, "5u9m-B.pdb"),
            Path(golden_data, "5u9m_updated.cif"),
            {
                "end": 248,
                "chain_id": "B",
                "pdb_id": "5u9m",
                "start": 1,
                "unp_end": 249,
                "coverage": 0.996,
                "unp_start": 2,
                "resolution": 2.35,
                "experimental_method": "X-ray diffraction",
                "tax_id": 559292,
            },
        ),
        (
            Path(golden_data, "5u9m-D.pdb"),
            Path(golden_data, "5u9m_updated.cif"),
            {
                "end": 248,
                "chain_id": "D",
                "pdb_id": "5u9m",
                "start": 1,
                "unp_end": 249,
                "coverage": 0.996,
                "unp_start": 2,
                "resolution": 2.35,
                "experimental_method": "X-ray diffraction",
                "tax_id": 559292,
            },
        ),
        (
            Path(golden_data, "1qup-A.pdb"),
            Path(golden_data, "1qup_updated.cif"),
            {
                "end": 222,
                "chain_id": "A",
                "pdb_id": "1qup",
                "start": 1,
                "unp_end": 223,
                "coverage": 0.892,
                "unp_start": 2,
                "resolution": 1.8,
                "experimental_method": "X-ray diffraction",
                "tax_id": 4932,
            },
        ),
        (
            Path(golden_data, "1qup-B.pdb"),
            Path(golden_data, "1qup_updated.cif"),
            {
                "end": 222,
                "chain_id": "B",
                "pdb_id": "1qup",
                "start": 1,
                "unp_end": 223,
                "coverage": 0.892,
                "unp_start": 2,
                "resolution": 1.8,
                "experimental_method": "X-ray diffraction",
                "tax_id": 4932,
            },
        ),
        (
            Path(golden_data, "1ej8-A.pdb"),
            Path(golden_data, "1ej8_updated.cif"),
            {
                "end": 140,
                "chain_id": "A",
                "pdb_id": "1ej8",
                "start": 1,
                "unp_end": 217,
                "coverage": 0.562,
                "unp_start": 78,
                "resolution": 1.55,
                "experimental_method": "X-ray diffraction",
                "tax_id": 4932,
            },
        ),
    ]


@pytest.fixture
def validate_api_hit_input():
    return [
        {
            "end": 249,
            "chain_id": "B",
            "pdb_id": "1jk9",
            "start": 1,
            "unp_end": 249,
            "coverage": 1,
            "unp_start": 1,
            "resolution": 2.9,
            "experimental_method": "X-ray diffraction",
            "tax_id": 4932,
        },
        {
            "end": 249,
            "chain_id": "D",
            "pdb_id": "1jk9",
            "start": 1,
            "unp_end": 249,
            "coverage": 1,
            "unp_start": 1,
            "resolution": 2.9,
            "experimental_method": "X-ray diffraction",
            "tax_id": 4932,
        },
        {
            "end": 248,
            "chain_id": "B",
            "pdb_id": "5u9m",
            "start": 1,
            "unp_end": 249,
            "coverage": 0.996,
            "unp_start": 2,
            "resolution": 2.35,
            "experimental_method": "X-ray diffraction",
            "tax_id": 559292,
        },
        {
            "end": 248,
            "chain_id": "D",
            "pdb_id": "5u9m",
            "start": 1,
            "unp_end": 249,
            "coverage": 0.996,
            "unp_start": 2,
            "resolution": 2.35,
            "experimental_method": "X-ray diffraction",
            "tax_id": 559292,
        },
        {
            "end": 222,
            "chain_id": "A",
            "pdb_id": "1qup",
            "start": 1,
            "unp_end": 223,
            "coverage": 0.892,
            "unp_start": 2,
            "resolution": 1.8,
            "experimental_method": "X-ray diffraction",
            "tax_id": 4932,
        },
        {
            "end": 222,
            "chain_id": "B",
            "pdb_id": "1qup",
            "start": 1,
            "unp_end": 223,
            "coverage": 0.892,
            "unp_start": 2,
            "resolution": 1.8,
            "experimental_method": "X-ray diffraction",
            "tax_id": 4932,
        },
        {
            "end": 140,
            "chain_id": "A",
            "pdb_id": "1ej8",
            "start": 1,
            "unp_end": 217,
            "coverage": 0.562,
            "unp_start": 78,
            "resolution": 1.55,
            "experimental_method": "X-ray diffraction",
            "tax_id": 4932,
        },
    ]


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
    pdb = selmodel_pdb(inp_pdb, 1)
    assert pdb.exists()
    pdb.unlink()


def test_validate_api_hit(
    mocker,
    fetch_pdb_files_output,
    pdb_hit_no_resolution,
    validate_api_hit_input,
):
    """Test validate_api_hit."""

    mock_fetch_pdb_files = mocker.patch("arctic3d.modules.pdb.fetch_pdb_files")
    mock_fetch_pdb_files.return_value = fetch_pdb_files_output

    validated_pdbs = validate_api_hit(
        fetch_list=validate_api_hit_input, uniprot_id="P40202"
    )
    assert validated_pdbs == fetch_pdb_files_output

    mock_fetch_pdb_files.return_value = [
        {
            "end": 140,
            "chain_id": "A",
            "pdb_id": "1ej8",
            "start": 1,
            "unp_end": 217,
            "coverage": 0.562,
            "unp_start": 78,
            "resolution": 1.55,
            "experimental_method": "X-ray diffraction",
            "tax_id": 4932,
        }
    ]
    validated_pdbs = validate_api_hit(
        fetch_list=validate_api_hit_input,
        uniprot_id="P40202",
        resolution_cutoff=1.6,
    )

    # Since `validate_api` will return the output of `fetch_pdb_files` (?)
    assert validated_pdbs == mock_fetch_pdb_files.return_value

    # Test the behaviour when the experimental method is NMR
    _validate_api_hit_input: dict[str, int | float | str | None] = (
        copy.deepcopy(validate_api_hit_input[0])
    )

    _validate_api_hit_input["experimental_method"] = "Solution NMR"
    _validate_api_hit_input["resolution"] = None
    _api_hit_input = [_validate_api_hit_input]

    mock_fetch_pdb_files.return_value = _api_hit_input

    validated_pdbs = validate_api_hit(
        fetch_list=_api_hit_input, uniprot_id="P20023"
    )

    # Since `validate_api` will return the output of `fetch_pdb_files` (?)
    assert validated_pdbs == mock_fetch_pdb_files.return_value


def test_get_best_pdb(
    mocker, fetch_pdb_files_output, inp_pdb_data, example_interfaces
):
    """Test get_best_pdb."""

    mock_make_request = mocker.patch("arctic3d.modules.pdb.make_request")
    mock_make_request.return_value = json.load(open(inp_pdb_data, "r"))

    mock_fetch_pdb_files = mocker.patch("arctic3d.modules.pdb.fetch_pdb_files")
    mock_fetch_pdb_files.return_value = fetch_pdb_files_output

    observed_pdb, observed_cif, filtered_interfaces = get_best_pdb(
        uniprot_id="P40202",
        interface_residues=example_interfaces,
    )  # type: ignore
    expected_pdb = Path("P40202-5u9m-B.pdb")
    expected_cif = Path("5u9m_updated.cif")
    exp_interfaces = {"P01024": [103, 104, 105]}

    mock_make_request.assert_called()
    mock_fetch_pdb_files.assert_called()

    assert observed_pdb is not None
    assert observed_cif is not None

    assert observed_pdb.name == expected_pdb.name
    assert observed_cif.name == expected_cif.name
    assert filtered_interfaces == exp_interfaces

    expected_pdb.unlink()

    # FIXME: What is this testing?
    orig_interfaces = {"P00441": [85, 137, 138]}
    observed_pdb, _, filtered_interfaces = get_best_pdb(
        "P40202", orig_interfaces, pdb_data=inp_pdb_data
    )  # type: ignore

    assert filtered_interfaces == orig_interfaces

    expected_pdb.unlink()


def test_fetch_pdb_files(mocker):

    mock_fetch_updated_cif = mocker.patch(
        "arctic3d.modules.pdb.fetch_updated_cif"
    )
    mock_fetch_updated_cif.return_value = Path("1ej8_updated.cif")

    mock_convert_cif_to_pdbs = mocker.patch(
        "arctic3d.modules.pdb.convert_cif_to_pdbs"
    )
    # `convert_cif_to_pdb` would produce a file, mock this
    converted_pdb = Path("1ej8-A.pdb")
    converted_pdb.touch()
    converted_cif = Path("1ej8_updated.cif")

    mock_convert_cif_to_pdbs.return_value = [converted_pdb]

    pdb_to_fetch = [
        {
            "end": 140,
            "chain_id": "A",
            "pdb_id": "1ej8",
            "start": 1,
            "unp_end": 217,
            "coverage": 0.562,
            "unp_start": 78,
            "resolution": 1.55,
            "experimental_method": "X-ray diffraction",
            "tax_id": 4932,
        }
    ]

    observed_output = fetch_pdb_files(
        pdb_to_fetch=pdb_to_fetch,
        uniprot_id="P40202",
    )
    expected_output = [(converted_pdb, converted_cif, pdb_to_fetch[0])]

    assert observed_output[0][0] == expected_output[0][0]
    assert observed_output[0][1] == expected_output[0][1]
    assert observed_output[0][2] == expected_output[0][2]

    converted_pdb.unlink()


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

    assert pdb_f is not None
    assert cif_f is not None
    assert top_hit is not None

    assert pdb_f.name == "4xoj-A-occ-tidy.pdb"
    assert cif_f.name == "4xoj_updated.cif"
    assert top_hit["pdb_id"] == "4xoj"
    assert top_hit["chain_id"] == "A"
    assert filtered_interfaces == {"P01024": [103, 104, 105]}

    pdb_f.unlink()
    cif_f.unlink()


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


def test_convert_cif_to_pdbs(inp_cif_3psg):
    """Test convert_cif_to_pdbs."""
    obs_out_pdb_fnames = convert_cif_to_pdbs(inp_cif_3psg, "3psg", "P00791")
    exp_out_pdb_fnames = [Path("3psg-A.pdb")]
    assert exp_out_pdb_fnames == obs_out_pdb_fnames
    # inspect the pdb file
    obs_pdb_lines = obs_out_pdb_fnames[0].read_text().splitlines()
    # checking first and last lines
    exp_pdb_lines = [
        "ATOM      1  N   LEU A  16      57.364  -9.595   2.554  1.00 21.58",
        "ATOM   2692  CB  ALA A 385      39.553 -10.495   1.923  1.00 22.44",
    ]
    assert obs_pdb_lines[0] == exp_pdb_lines[0]
    assert obs_pdb_lines[-1] == exp_pdb_lines[-1]

    obs_out_pdb_fnames[0].unlink()
