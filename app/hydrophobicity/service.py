from .constants import HYDROPHOBICITY_SCALE

#|------------------------------------------------------------------------------|#


class ProteinHydrophobicity():
    def __init__(self, sequence):
        self.sequence = sequence

    #|------------------------------------------------------------------------------|#

    def hydrophobicity(self, scale = "Kyte-Doolittle") -> float:

        hydrophobicity_scale = HYDROPHOBICITY_SCALE[scale]
        hydrophobicity = 0.0

        for aa in self.sequence:
            hydrophobicity += hydrophobicity_scale[aa]

        return hydrophobicity / len(self.sequence)
        