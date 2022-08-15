"""
Created by Daniel Fu (Reif Lab, Duke University) at 6/18/2020

Name        : solve_crossovers.py
Project     : cadaxisdna
Description : Solves the optimum placement of all crossovers
                See, https://www.geeksforgeeks.org/place-k-elements-such-that-minimum-distance-is-maximized/
                # Python 3 program to find largest minimum
                # distance among k points.
Interpreter : Python 3.7.4
"""

# from ..routing.helper import log
import math


def isFeasible(mid, arr, n, k):
    """
    Returns true if it is possible to arrange
    k elements of arr[0..n-1] with minimum
    distance given as mid.
    :param mid: Minimum accepted distance
    :param arr: Elements of the array
    :param n: Length of array
    :param k: Number of elements to place
    :return: True is elements fit, False otherwise
    """
    # Place first element at arr[0] position
    pos = arr[0]

    # Initialize count of elements placed.
    elements = 1

    # Try placing k elements with minimum
    # distance mid.
    for i in range(1, n, 1):
        if arr[i] - pos >= mid:

            # Place next element if its distance
            # from the previously placed element
            # is greater than current mid
            pos = arr[i]
            elements += 1

            # Return if all elements are placed
            # successfully
            if elements == k:
                return True
    return False


def get_spacing(arr, n, spacing, k):
    """
    Calculates the array given a spacing
    :param arr: input array
    :param n: length of array
    :param spacing: desired distance between elements
    :param k: number of elements
    :return:
    """
    arr_out = []

    pos = arr[0]
    arr_out.append(pos)
    for i in range(1, n, 1):
        if arr[i] - pos >= spacing:
            pos = arr[i]
            arr_out.append(pos)
            if len(arr_out) == k:
                return arr_out
    return arr_out


def largestMinDist(arr, k):
    """
    # Returns largest minimum distance for k elements
    # in arr[0..n-1]. If elements can't be placed,
    # returns -1.

    All elements of array must be positive
    :param arr: position of elements
    :param k: number of elements
    :return:
    """
    # Sort the positions
    arr.sort(reverse=False)

    # Catch edge case
    if k == 1:
        pos = math.ceil(len(arr) / 2)
        res = min(arr[-1] - pos, pos - arr[0])
        return [arr[pos]], res, True

    # Initialize result.
    res = -1

    # Determine the size of the array
    n = len(arr)

    # Consider the maximum possible distance
    left = 0
    right = arr[n - 1]

    # Do binary search for largest
    # minimum distance
    while right - left > 0.001:
        mid = (left + right) / 2

        # If it is possible to place k elements
        # with minimum distance mid, search for
        # higher distance.
        if isFeasible(mid, arr, n, k):
            # Change value of variable max to mid iff
            # all elements can be successfully placed
            res = max(res, mid)
            left = mid

        # If not possible to place k elements,
        # search for lower distance
        else:
            res = mid
            right = mid
    res = round(res, 1)
    arr_out = get_spacing(arr, n, res, k)
    if len(arr_out) == k:
        valid = True
    else:
        valid = False
    return arr_out, res, valid


def iterate(arr_in, mindist):
    """
    Finds largest array from input elements (arr_in) that satisfy a minimum distance spacing (mindist) by
    iterating through resulting number of elements
    :param arr_in:
    :param mindist:
    :return:
    """
    n = len(arr_in)
    for k in range(1, n):
        arr, spacing, valid = largestMinDist(arr_in, k)
        # log.system("Largest set is", arr_out, spacing, valid)
        if spacing > mindist and valid:
            arr_out = arr
            spacing_out = spacing
        else:
            pass
    try:
        # Edge case, for all n, spacing < mindist or invalid, so arr_out never set
        return arr_out, spacing_out
    except UnboundLocalError:
        return [], -1


# Driver code
if __name__ == '__main__':
    import random
    # ### This code tests iterate
    # # Array with the integer positions of each crossover, assumng the first crossover is at position 0
    # arr = [random.randint(0, 400) for x in range(int(400/10.5))]
    # print(arr)
    # # Desired minimum distance
    # mindist = 21
    # print(iterate(arr, mindist))

    ### This code tests largestMinDist
    # Array with the integer positions of each crossover, assumng the first crossover is at position 0
    # arr = [random.randint(0, 400) for x in range(int(400 / 10.5))]
    arr = [26, 36, 56, 76, 96, 116]
    print(arr)
    # Number of crossovers to place
    k = 4
    print(largestMinDist(arr, k))


# This code is contributed by
# Sanjit_prasad
