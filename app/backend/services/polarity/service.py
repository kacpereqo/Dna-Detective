from .constants import POLARITY_SCALE

# |------------------------------------------------------------------------------|#


class ProteinPolarity():
    def __init__(self, sequence):
        self.sequence = sequence

    # |------------------------------------------------------------------------------|#

    def get_avg_polarity(self) -> float:

        polarity = 0.0

        for aa in self.sequence:
            polarity += POLARITY_SCALE[aa]

        return polarity / len(self.sequence)

    # |------------------------------------------------------------------------------|#

    def get_polarity(self, window: int = 3):
        data = {}
        start = int(window / 2)

        for aa in range(start, len(self.sequence) - start):
            polarity = round(sum([POLARITY_SCALE[self.sequence[aa + i]]
                                  for i in range(-start, start + 1)]) / window, 1)
            data[aa] = polarity
        return data
