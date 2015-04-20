#!/usr/bin/env python
if __name__ == '__main__':
    import sys
    import os

    #python magic to add the current directory to the pythonpath
    sys.path.append(os.getcwd())

    # ==========================================================
    # Setup django environment 
    # ==========================================================
    if not os.environ.has_key('DJANGO_SETTINGS_MODULE'):
        os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'
    # ==========================================================
    # Done setting up django environment
    # ==========================================================


import math
from pgd.pgd_core.models import *
import json

def dump_protein(code, preserve_precision):
    """
    Dumps a protein to JSON so that it may be imported easily into
    another program for comparison
    """
    protein = Protein.objects.get(code=code)

    precision = 3
    convert_str = '%%.%if' % (precision)

    count = 0

    # A list of values that should not be printed out
    FIELDS = {'aa':-1,
            'ss':-1,
            'a1':1,
            'a2':1,
            'a3':1,
            'a4':1,
            'a5':1,
            'a6':1,
            'a7':1,
            'L1':3,
            'L2':3,
            'L3':3,
            'L4':3,
            'L5':3,
            'phi':2,
            'psi':2,
            'ome':2,
            'chi':2,
            'bm':3,
            'bs':3,
            'bg':3, 
            'h_bond_energy':3,
            'zeta':2}

    _json = {'code':code,
            'resolution':protein.resolution,
            'threshold':protein.threshold,
            'rfactor':protein.rfactor
            }

    json_residues = {}
    for r in protein.residues.all():
        json_residue = {}
        for field, precision in FIELDS.items():
            if preserve_precision or precision == -1:
                json_residue[field] = r.__dict__[field]

            else:
                precision_str = '%%.%if' % precision
                json_residue[field] = precision_str % r.__dict__[field]

        json_residues[int(r.oldID)] = json_residue

    _json['residues'] = json_residues

    return json.dumps(_json)


if __name__ == '__main__':
    import sys

    try:
        code = sys.argv[1]
        for code in sys.argv[1:]:
            print dump_protein(code, True)
    except IndexError:
        print 'Usage: dump_protein.py code [code2 [code3 [...]]] '