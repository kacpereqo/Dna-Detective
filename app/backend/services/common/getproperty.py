def get_property(scale, sequence):
    data = {}
    window = 3
    start = int(window / 2)

    for aa in range(start, len(sequence) - start):
        polarity = round(sum([scale[sequence[aa + i]]
                              for i in range(-start, start + 1)]) / window, 1)
        data[aa] = polarity
    return data
