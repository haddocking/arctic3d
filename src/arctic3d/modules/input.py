class Input:
    def __init__(self, input_arg):
        self.arg = input_arg

    def is_fasta(self):
        return self.arg.endswith(".fasta")

    def is_pdb(self):
        return self.arg.endswith(".pdb")

    def is_uniprot(self):
        return len(self.arg.split(".")) == 1
