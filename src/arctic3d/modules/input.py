def is_fasta(input_fasta):
    if input_fasta:
        return input_fasta.endswith(".fasta")
    else:
        return False


def is_pdb(input_pdb):
    if input_pdb:
        return input_pdb.endswith(".pdb")
    else:
        return False


def is_uniprot(input_uniprot):
    if input_uniprot:
        return len(input_uniprot.split(".")) == 1
    else:
        return False
