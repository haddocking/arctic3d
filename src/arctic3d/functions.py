import time
from typing import Union

import requests


def make_request(
    url: str,
    data: Union[
        dict[str, list[dict[str, Union[str, int, float, None]]]], None
    ],
) -> Union[dict[str, list[dict[str, Union[str, int, float, None]]]], None]:
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
