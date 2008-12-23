""" *******************************************************************************
    Prints out the query results in a HTML format or in a dump file
    in an HTML format the function will paginate the data
    numResidues = total # of residues in query
    Result_Set = set of results to be printed
    Per_Page = number of results to show per page
    toFile = file to dump results to
    **************************************************************************  """

import math
from pgd_search.models import *
from constants import AA_CHOICES

def dumpSearch(search, writer):

    length   = search.segmentLength
    querySet = search.querySet()

    count = 0

    # A list of values that should not be printed out
    FIELDS = ['aa','a1','a2','a3','a4','a5','a6','a7','L1','L2','L3','L4','L5','ss','phi', 'psi', 'ome', 'chi', 'bm', 'bs', 'bg', 'h_bond_energy', 'zeta']
    FIELD_LABEL_REPLACEMENTS = {'h_bond_energy':'H Bond', 'aa':'AA'}
    FIELD_VALUE_REPLACEMENTS = {'aa':AA_CHOICES}

    total = querySet.count()

    # General query result information
    writer.write("Match\tCode\tResidue\tID")

    # Field names
    for field in FIELDS:
        writer.write('\t')
        if field in FIELD_LABEL_REPLACEMENTS:
            writer.write(FIELD_LABEL_REPLACEMENTS[field])
        else:
            writer.write(field)
    writer.write('\n')

    #calculate list of iValues
    iValues = []
    start = 0 - (length-1) / 2
    stop  = int(math.ceil((length-1) / 2.0))+1
    for i in range(start, stop):
        if i < 0:
            iValues.append((i,'(i - %i)' % math.fabs(i)))
        elif i == 0:
            iValues.append((i,'(i)'))
        else:
            iValues.append((i,'(i+%i)'% i))

    #calculate iIndex
    iIndex = int(math.ceil(length/2.0)-1)

    # Go Time. Print out the data.
    for segment in querySet:
        count += 1
        first = True
        for offset, string in iValues:

            #segment count
            if first:
                writer.write(count)
                first=False
            writer.write('\t')

            #protein code
            writer.write(segment.protein_id)
            writer.write('\t')

            #index string repr
            writer.write(string)
            writer.write('\t')

            #residue index
            writer.write(segment.__dict__['r%i_chainIndex' % (iIndex+offset) ])

            #field values
            for field in FIELDS:
                writer.write('\t')
                # replace field with display value if needed
                if field in FIELD_VALUE_REPLACEMENTS:
                    code = segment.__dict__['r%i_%s' % (iIndex+offset, field)]
                    if code:
                        for k,v in FIELD_VALUE_REPLACEMENTS[field]:
                            if k == code:
                                writer.write(v)
                # just write value
                else:
                    writer.write(segment.__dict__['r%i_%s' % (iIndex+offset, field)])

            #end residue with a newline
            writer.write('\n')