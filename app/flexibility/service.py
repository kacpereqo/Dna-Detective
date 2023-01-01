from .constants import FLEXIBILITY

class ProteinFlexibility():
    def __init__(self, protein):
        self.protein = protein

    #|------------------------------------------------------------------------------|#

    def flexibility(self, window_size = 9) -> float:
        """Returns the flexibility of the protein"""
        weights = (0.25, 0.4375, 0.625, 0.8125, 1)
        scores = []

        for i in range(self.length - window_size):
            subsequence = self.sequence[i : i + window_size]
            score = 0.0

            for j in range(window_size // 2):
                front = subsequence[j]
                back = subsequence[window_size - j - 1]
                score += (FLEXIBILITY[front] + FLEXIBILITY[back]) * weights[j]

            middle = subsequence[window_size // 2 + 1]
            score += FLEXIBILITY[middle]

            scores.append(score / 5.25)

        return scores
