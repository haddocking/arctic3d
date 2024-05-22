from typing import Any

from clustering import interface_clustering

interface_dict: dict[str, Any] = {
    "int_1": [],
    "int_2": [],
    "int_3": [],
    "int_4": [],
}

interface_clustering(
    interface_dict, "../../../tests/golden_data/interface_matrix.txt"
)
