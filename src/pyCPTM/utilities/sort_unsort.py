"""Return a sorted array and the indices to unsort the array, i.e. reverse the sorting of the original array."""

import numpy as np


def sort_unsort(a):
    """
    Return a sorted array and the indices to unsort the array, i.e. reverse the sorting of the original array.

    Parameters
    ----------
    a : ndarray
        1D array to be sorted.

    Returns
    -------
    sortedArray : ndarray
        Sorted array
    unsortIndex : ndarray
        Array of indices to unsort/reverse the sorting. Use as sortedArray[unsortIndex] to get 'a'.

    """
    if not isinstance(a, np.ndarray):
        print("Input array must be np.ndarray")
    else:
        # Sort array and get sorting indices
        sortedArray = np.sort(a)
        sortedArrayId = np.argsort(a)
        # Indices to "unsort" the sorted array
        unsortIndex = np.argsort(sortedArrayId)
        unsortedArray = sortedArray[unsortIndex]
        if all(unsortedArray == a):
            return sortedArray, unsortIndex
        else:
            print("Sorting failed.")
