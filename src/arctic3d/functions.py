import time

import requests


def make_request(
    url: str, data: dict[str, list[dict[str, str | int | float | None]]] | None
) -> dict[str, list[dict[str, str | int | float | None]]] | None:
    """Helper function to make the requests."""
    for _ in range(3):
        response = requests.get(url)
        if response.status_code != 404:
            data = response.json()
            break
    if response.status_code == 404:
        data = None
    time.sleep(0.1)
    return data
