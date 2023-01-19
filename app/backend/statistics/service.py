class ProteinStatistics():
    def __init__(self, protein):
        self.protein = protein

    #|------------------------------------------------------------------------------|#

    def perctage_of(self, amino_acid: str) -> float:
        """Returns the percentage of a given amino acid in the protein"""
        return round((self.protein.count(amino_acid) / len(self.protein)) * 100, 2)

    #|------------------------------------------------------------------------------|#

    def percentage_of_all(self) -> dict:
        """Returns the percentage of all amino acids in the protein"""
        return {amino_acid: self.perctage_of(amino_acid) for amino_acid in self.protein}

#|------------------------------------------------------------------------------|#

    def count_of(self, amino_acid: str) -> int:
        """Returns the count of a given amino acid in the protein"""
        return self.protein.count(amino_acid)

#|------------------------------------------------------------------------------|#

    def count_of_all(self) -> dict:
        """Returns the count of all amino acids in the protein"""
        return {amino_acid: self.count_of(amino_acid) for amino_acid in self.protein}

#|------------------------------------------------------------------------------|#

    def percentage_of_charged(self) -> dict:
        """Returns the percentage of charged amino acids in the protein"""
        charged = ["R", "H", "K", "D", "E"]
        return {amino_acid: self.perctage_of(amino_acid) for amino_acid in charged}

#|------------------------------------------------------------------------------|#

    def count_of_charged(self) -> dict:
        """Returns the count of charged amino acids in the protein"""
        charged = ["R", "H", "K", "D", "E"]
        return {amino_acid: self.count_of(amino_acid) for amino_acid in charged}

#|------------------------------------------------------------------------------|#

    