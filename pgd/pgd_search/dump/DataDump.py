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
from pgd_constants import AA_CHOICES
from pgd_search.models import searchSettings
from django.core.paginator import Paginator

def dumpSearch(search, writer):

    length   = search.segmentLength
    querySet = search.querySet()

    count = 0

    # A list of values that should not be printed out
    FIELDS = ['aa','a1','a2','a3','a4','a5','a6','a7','L1','L2','L3','L4','L5','ss','phi', 'psi', 'ome', 'chi', 'bm', 'bs', 'bg', 'h_bond_energy', 'zeta']
    FIELD_LABEL_REPLACEMENTS = {
        'h_bond_energy':'H Bond', 
        'aa':'AA',
        'L1':u'C(-1)N',
        'L2':u'N-CA',
        'L3':u'CA-CB',
        'L4':u'CA-C',
        'L5':'C-O',
        'a1':u'C(-1)-N-CA',
        'a2':u'N-CA-CB',
        'a3':u'N-CA-C',
        'a4':u'CB-CA-C',
        'a5':u'CA-C-O',
        'a6':u'CA-C-N(+1)',
        'a7':u'O-C-N(+1)'
        }
    FIELD_VALUE_REPLACEMENTS = {'aa':AA_CHOICES}

    total = querySet.count()

    # General query result information
    writer.write("Match\tCode\tResidue\tID\tChain ID")

    # Field names
    for field in FIELDS:
        writer.write('\t')
        if field in FIELD_LABEL_REPLACEMENTS:
            writer.write(FIELD_LABEL_REPLACEMENTS[field])
        else:
            writer.write(field)
    writer.write('\n')

    #calculate list of iValues
    iValues = [
        (
            i, #index int
            ('(i%+i)' % i) if i else '(i)', # string representation
        ) for i in range(
            0 - (length-1)/2, #start
            int(math.ceil((length-1) / 2.0))+1, #stop
        )
    ] 

    #calculate iIndex
    iIndex = int(math.ceil(searchSettings.segmentSize/2.0)-1)

    # Go Time. Print out the data.
    chunks = Paginator(querySet,500)
    for pageNo in chunks.page_range:
        for segment in chunks.page(pageNo).object_list:
            count += 1
            first = True
            for offset, string in iValues:

                #segment count
                if first:
                    writer.write(count)
                    first=False
                else:
                    writer.write(' ')
                writer.write('\t')

                #protein code
                writer.write(segment.protein_id)
                writer.write('\t')

                #index string repr
                writer.write(string)
                writer.write('\t')

                #residue index
                writer.write(segment.__dict__['r%i_oldID' % (iIndex+offset) ])
                writer.write('\t')

                #chain Index
                writer.write(segment.chainID)

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
