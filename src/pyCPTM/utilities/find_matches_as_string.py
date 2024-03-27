import numpy as np


def find_matches_as_string(arrayA, arrayB):
    matches = set(arrayA) & set(arrayB)
    matches = np.array(list(matches), dtype=np.dtype("U"))
    matches.sort()
    return matches
