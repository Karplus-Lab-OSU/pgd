import math


def residue_indexes(length):
    """ returns a list of i-indexes for the given length """
    return range(0 - (length-1) / 2, int(math.ceil((length-1) / 2.0))+1)
