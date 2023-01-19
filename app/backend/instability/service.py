from Bio.SeqUtils.ProtParam import ProteinAnalysis

class ProteinInstability():
    def __init__(self, protein):
        self.protein = protein

    #|------------------------------------------------------------------------------|#

    def instability(self) -> float:
        """Returns the flexibility of the protein"""
        sequence = ProteinAnalysis(self.protein)
        return sequence.instability_index()
        
