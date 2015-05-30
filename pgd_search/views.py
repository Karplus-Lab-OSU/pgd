import math

from django.conf import settings


"""
Properties for residues as indexes
"""
RESIDUE_INDEX_START = 0 - (settings.SEGMENT_SIZE-1) / 2
RESIDUE_INDEX_STOP = int(math.ceil((settings.SEGMENT_SIZE-1) / 2.0))+1
RESIDUE_INDEXES = range(RESIDUE_INDEX_START, RESIDUE_INDEX_STOP)


def settings_processor(request):
    """
    settings_processor adds settings required by most pages
    """

    return {
        'PGD_VERSION': settings.PGD_VERSION,
        'ROOT': settings.SITE_ROOT,
        'DATA_VERSION': settings.DATA_VERSION,
        'GOOGLE_ID': settings.GOOGLE_ID
    }
