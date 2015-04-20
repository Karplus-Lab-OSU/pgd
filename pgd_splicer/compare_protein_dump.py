#!/usr/bin/env python
import json

def compare_dump(file_new, file_old):

    #load files and convert them to dicts
    _file = None
    try:
        _file = open(file_new,'r')
        new = json.loads(_file.read())

        _file = open(file_old,'r')
        old = json.loads(_file.read())

    finally:
        if _file:
            _file.close()


    print 'comparing proteins.  fields are always listed: New, Old'

    # compare protein properties
    print 'Protein Properties:'
    properties = ['code','threshold','resolution','rfactor']
    for prop in properties:
        if new[prop] != old[prop]:
            print '   %s: %s != %s' % (prop, new[prop], old[prop])


    # Compare the id's checking to make sure the lists match.  If a residue is
    # missing from either list we'll remove it.  when we iterate the residues
    # to compare properties we'll have only matching residues
    print 'Missing Residues:'
    missing = False
    for index in new['residues'].keys():
        if not index in old['residues']:
            missing = True
            print 'Missing Residue (old): ', index
            del new['residues'][index]

    for index in old['residues'].keys():
        if not index in new['residues']:
            missing = True
            print 'Missing Residue (new): ', index
            del old['residues'][index]

    if not missing:
        print '   None'


    # Compare properties of all residues
    print 'Residue Properties:'
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
            #'h_bond_energy':3,
            'zeta':2}

    keys = sorted([int(i) for i in new['residues']])

    for i in keys:
        key = str(i)
        rn = new['residues'][key]
        ro = old['residues'][key]
        for field, precision in FIELDS.items():
            if precision == -1:
                if field == 'ss':
                    if rn[field] != ro[field] and not(rn[field]=='-' and ro[field] == ' '):
                        print '    [%s %s] %s: %s != %s' % (i, rn['aa'], field, rn[field], ro[field])
                elif rn[field] != ro[field]:
                    print '    [%s %s] %s: %s !=%s' % (i, rn['aa'], field, rn[field], ro[field])

            else:
                precision_str = '%%.%if' % precision
                formated = precision_str % rn[field]
                compare = not formated == ro[field] \
                          and not (rn[field] == 999.9 and float(ro[field]) == 0.0) \
                          and not ((float(formated)+0.001).__str__() == ro[field] or (float(formated)-0.001).__str__() == ro[field])
                if compare:
                    output = '    [%s %s] %s: %s ('+precision_str+') != %s'
                    print output % (i, rn['aa'], field, rn[field], rn[field], ro[field])


def compare_fields(new, old, precision):
    """
    Function for comparing fields.  this is its own function
    in case we need to use some complicated logic to deal with the different
    precision between the numbers
    """
    precision_str = '%%.%if' % precision
    return precision_str % new == old


if __name__ == '__main__':
    import sys

    try:
        new = sys.argv[1]
        old = sys.argv[2]
        compare_dump(new, old)

    except IndexError:
        print 'Usage: compare_protein_dump.py new_file old_file'
        print ''
        print '    new_file - file from new splicer'
        print '    old_file - file from old splicer'