"""Tool to create a benchmark set for the ARCTIC-3D tool."""
import argparse
import logging
import os
from typing import Optional

import pandas as pd  # type: ignore
import requests

logging.basicConfig(
    format=" %(funcName)s:L%(lineno)d %(levelname)s - %(message)s",
    level=logging.INFO,
    datefmt="%m/%d/%Y %I:%M:%S %p",
)


def load_bm5_file(bm5_file: str) -> list[tuple[str, str]]:
    """Load the BM5 Excel File."""
    logging.info("Loading BM5 file.")
    bm5_df = pd.read_excel(bm5_file, header=2).dropna()

    receptor_ids = bm5_df["PDB ID 1"].values
    ligand_ids = bm5_df["PDB ID 2"].values

    data = []
    for receptor, ligand in zip(receptor_ids, ligand_ids):
        data.append((receptor, ligand))

    return data


def filter_bm5(bm5: list[tuple[str, str]]) -> list[tuple[str, str]]:
    """Filter out multichain pdbs."""
    logging.info("Filtering out multichain pdbs.")
    filtered_bm5 = []
    for element in bm5:
        receptor, ligand = element
        receptor_chains = len(receptor.split("_")[-1])
        ligand_chains = len(ligand.split("_")[-1])
        if receptor_chains != 1 or ligand_chains != 1:
            continue
        filtered_bm5.append(element)
    return filtered_bm5


def identify_uniprotid(pdb_id: str, target_chain: str) -> Optional[str]:
    """Identify the Uniprot ID for a given PDB ID."""
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
    response = requests.get(url)

    if response.status_code != 200:
        logging.warning(f"No Uniprot annotation for {pdb_id}_{target_chain}.")
        return None

    uniprot_data = response.json()[pdb_id.lower()]["UniProt"]

    for uniprot_id in uniprot_data:
        for entry in uniprot_data[uniprot_id]["mappings"]:
            chain = entry["chain_id"]
            if chain == target_chain:
                return uniprot_id
    return None


def parse_bm(bm: list[tuple[str, str]]) -> list[tuple[str, str, str, str]]:
    """Parse the benchmark and return a list with the Uniprot IDs."""
    parsed_bm: list = []
    for receptor, ligand in bm:
        receptor_pdb, receptor_chain = receptor.split("_")
        ligand_pdb, ligand_chain = ligand.split("_")
        uniprot_receptor = identify_uniprotid(receptor_pdb, receptor_chain)
        uniprot_ligand = identify_uniprotid(ligand_pdb, ligand_chain)

        if uniprot_receptor and uniprot_ligand:
            parsed_bm.append((receptor, uniprot_receptor, ligand, uniprot_ligand))

    return parsed_bm


def write_output_file(
    parsed_bm: list[tuple[str, str, str, str]], output_file: str
) -> None:
    """Write the output file."""
    logging.info(f"Saving output to {output_file}.")
    with open(output_file, "w") as f:
        f.write("receptor,uniprot_receptor,ligand,uniprot_ligand" + os.linesep)
        for receptor, uniprot_receptor, ligand, uniprot_ligand in parsed_bm:
            f.write(
                f"{receptor},{uniprot_receptor},{ligand},{uniprot_ligand}" + os.linesep
            )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument("output_file")
    args = parser.parse_args()

    bm5_input_file = args.input_file
    output_file = args.output_file

    bm_list = load_bm5_file(bm5_input_file)
    filtered_bm_list = filter_bm5(bm5=bm_list)
    parsed_bm = parse_bm(filtered_bm_list)
    write_output_file(parsed_bm, output_file)


if __name__ == "__main__":
    main()
