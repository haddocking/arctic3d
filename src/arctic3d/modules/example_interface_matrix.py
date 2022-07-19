from interface_matrix import interface_matrix

interface_dict = {
    "int_1": [16, 18, 19],
    "int_2": [21, 22, 24, 26, 27, 29, 30, 31, 36, 38, 39, 40],
    "int_3": [16, 17, 18],
    "int_4": [20, 24],
    "int_5": [1, 2, 3],
    "int_6": [1, 2],
}

interface_matrix(interface_dict, "../../../tests/golden_data/1rypB_r_b.pdb")
