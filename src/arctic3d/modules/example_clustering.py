from clustering import interface_clustering

interface_dict = {
    "int_1": [],
    "int_2": [],
    "int_3": [],
    "int_4": [],
}

interface_clustering(
    interface_dict, "../../../tests/golden_data/interface_matrix.txt"
)
