import math

from pgd_search.models import searchSettings
from pgd_splicer.models import pdb_select_settings
import settings

"""
Properties for residues as indexes
"""
RESIDUE_INDEX_START = 0 - (searchSettings.segmentSize-1) / 2
RESIDUE_INDEX_STOP  = int(math.ceil((searchSettings.segmentSize-1) / 2.0))+1
RESIDUE_INDEXES = range(RESIDUE_INDEX_START,RESIDUE_INDEX_STOP)


def settings_processor(request):
    """
    settings_processor adds settings required by most pages
    """

    return {
        'PGD_VERSION':settings.PGD_VERSION,
        'MEDIA':settings.MEDIA_URL,
        'ROOT':settings.SITE_ROOT,
        'DATA_VERSION':pdb_select_settings.DATA_VERSION
    }