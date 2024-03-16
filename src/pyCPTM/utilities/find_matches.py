# -*- coding: utf-8 -*-
import numpy as np


def find_matches(arrayA, arrayB):
    matches = set(arrayA) & set(arrayB)
    matches = np.array(list(matches))
    matches.sort()
    return matches
