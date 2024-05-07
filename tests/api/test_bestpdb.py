from arctic3d.modules.pdb import BESTPDB_URL
import requests
import http
from pathlib import Path
import json
import pytest

TARGET_UNIPROT = "P07550"

EXPECTED_RESPONSE = Path(Path(__file__).parent, "bestpdb.json")


@pytest.mark.sanity
def test_bestpdb(compare_responses: callable):
    # Make a request to the endpoint and check if the response is 200
    response = requests.get(f"{BESTPDB_URL}/{TARGET_UNIPROT}")
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
    with open(EXPECTED_RESPONSE, "r") as f:
        expected_response = json.load(f)

    compare_responses(
        response.json(), expected_response
    ), "Observed response is missing keys, the API might have changed"
