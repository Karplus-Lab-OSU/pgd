import math
from pgd_search.models import searchSettings

"""
Properties for residues as indexes
"""
RESIDUE_INDEX_START = 0 - (searchSettings.segmentSize-1) / 2
RESIDUE_INDEX_STOP  = int(math.ceil((searchSettings.segmentSize-1) / 2.0))+1
RESIDUE_INDEXES = range(RESIDUE_INDEX_START,RESIDUE_INDEX_STOP)
