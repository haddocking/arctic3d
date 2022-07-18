from interface import get_interface_residues

# first example : no filters
filters_dict = {
    "outHM": False,
    "outHD": False,
    "cpc": False,
    "cbi": False,
    "outIg": False,
    "t2m": False,
}

uniprot_ids = ["P49916", "P24941", "P00690"]
for uniprot_id in uniprot_ids:
    output_data = get_interface_residues(uniprot_id, filters_dict)
    print(f"interface data for {uniprot_id} : {output_data}")
