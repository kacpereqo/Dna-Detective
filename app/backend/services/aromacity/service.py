from Bio.SeqUtils.ProtParam import ProteinAnalysis

class ProteinAromacity():
    def __init__(self, protein):
        self.protein = protein

    #|------------------------------------------------------------------------------|#

    def aromacity(self) -> float:
        """Returns the flexibility of the protein"""
        sequence = ProteinAnalysis(self.protein)
        return sequence.aromaticity()
        
