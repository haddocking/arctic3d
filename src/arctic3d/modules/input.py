class Input:
    def __init__(self, input_arg):
        self.arg = input_arg

    def is_fasta(self):
        return self.arg.endswith(".fasta")
