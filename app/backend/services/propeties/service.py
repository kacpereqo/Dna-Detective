from .constants import AMINE_ACIDS_MASS


class SequncePropeties():
    def __init__(self, sequence: str) -> None:
        self.sequence = sequence

    def get_weight(self, window: int = 3):
        weights = []
        start = int(window / 2)

        for aa in range(start, len(self.sequence) - start):
            weight = round(sum([AMINE_ACIDS_MASS[self.sequence[aa + i]]
                                for i in range(-start, start + 1)]) / window, 1)
            weights.append(weight)
        return weights
