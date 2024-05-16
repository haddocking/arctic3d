import http
import json
from pathlib import Path
from typing import Callable

import pytest
import requests

from arctic3d.cli_localise import UNIPROT_API_URL
from arctic3d.modules.interface import LIGAND_URL
from arctic3d.modules.pdb import BESTPDB_URL, PDBE_URL

from . import golden_data

TARGET_UNIPROT = "P07550"
TARGET_PDB = "1crn"


@pytest.mark.sanity
def test_bestpdb(compare_responses: Callable):
    # Make a request to the endpoint and check if the response is 200
    response = requests.get(f"{BESTPDB_URL}/{TARGET_UNIPROT}", timeout=120)
    assert (
        response.status_code == http.HTTPStatus.OK
    ), f"Endpoint {BESTPDB_URL} not reachable"

    # Check if the response is a JSON
    assert (
        response.headers["Content-Type"] == "application/json"
    ), "Response is not JSON"

    # Check if the response is not empty
    assert response.json(), "Response is empty"

    # Check if the response has the same fields as the expected response
    # Load the expected response
    with open(Path(golden_data, "bestpdb.json"), "r") as f:
        expected_response = json.load(f)

    assert compare_responses(
        response.json(), expected_response
    ), "Observed response is missing keys, the API might have changed"


@pytest.mark.sanity
def test_pdbe():
    response = requests.get(f"{PDBE_URL}/pdb{TARGET_PDB}.ent")

    assert (
        response.status_code == http.HTTPStatus.OK
    ), f"Endpoint {PDBE_URL} not reachable"

    # Check if the response is a text file
    assert (
        response.headers["Content-Type"] == "text/plain; charset=UTF-8"
    ), "Response is not a PDB file"

    # Check if the response is not empty
    assert response.text, "Response is empty"

    response = requests.get(f"{PDBE_URL}/{TARGET_PDB}_updated.cif")

    assert (
        response.status_code == http.HTTPStatus.OK
    ), f"Endpoint {PDBE_URL} not reachable"

    # Check if the response is a text file
    assert (
        response.headers["Content-Type"] == "text/plain; charset=UTF-8"
    ), "Response is not a CIF file"

    # Check if the response is not empty
    assert response.text, "Response is empty"


@pytest.mark.sanity
def test_uniprot(compare_responses):
    response = requests.get(f"{UNIPROT_API_URL}/{TARGET_UNIPROT}")

    assert (
        response.status_code == http.HTTPStatus.OK
    ), f"Endpoint {UNIPROT_API_URL} not reachable"

    # Check if the response is a JSON
    assert (
        response.headers["Content-Type"] == "application/json"
    ), "Response is not JSON"

    # Check if the response is not empty
    assert response.json(), "Response is empty"

    # Check if the response has the same fields as the expected response
    # Load the expected response
    with open(Path(golden_data, "uniprot.json"), "r") as f:
        expected_response = json.load(f)

    assert compare_responses(
        response.json(), expected_response
    ), "Observed response is missing keys, the API might have changed"


@pytest.mark.sanity
def test_ligand():
    response = requests.get(f"{LIGAND_URL}/{TARGET_UNIPROT}")

    assert (
        response.status_code == http.HTTPStatus.OK
    ), f"Endpoint {LIGAND_URL} not reachable"

    # Check if the response is a JSON
    assert (
        response.headers["Content-Type"] == "application/json"
    ), "Response is not JSON"

    # Check if the response is not empty
    assert response.json(), "Response is empty"

    # Check if the response has the same fields as the expected response
    # Load the expected response
    with open(Path(golden_data, "ligand.json"), "r") as f:
        expected_response = json.load(f)

    assert (
        response.json() == expected_response
    ), "Observed response is different from the expected response"
