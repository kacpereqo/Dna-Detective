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