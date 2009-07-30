from django.shortcuts import render_to_response
from django.template import RequestContext
from django.conf import settings
import math

from pgd_constants import AA_CHOICES
from pgd_search.models import Search, Segment, iIndex
from pgd_search.plot.ConfDistFuncs import getCircularStats

ANGLES_BASE = ('ome', 'phi', 'psi', 'chi', 'zeta')
angles = ['r%i_%s' %(iIndex, angle) for angle in ANGLES_BASE]

"""
display statistics about the search
"""
def searchStatistics(request):

    #store globals locally for speed optimization

    stat_attributes_base = [('L1',u'C\u207B\u00B9N'),
                        ('L2',u'NC\u1D45'),
                        ('L3',u'C\u1D45C\u1D5D'),
                        ('L4',u'C\u1D45C'),
                        ('L5','CO'),
                        ('a1',u'C\u207B\u00B9NC\u1D45'),
                        ('a2',u'NC\u1D45C\u1D5D'),
                        ('a3',u'NC\u1D45C'),
                        ('a4',u'C\u1D5DC\u1D45C'),
                        ('a5',u'C\u1D45CO'),
                        ('a6',u'C\u1D45CN\u207A\u00B9'),
                        ('a7',u'OCN\u207A\u00B9'),
                        ('ome',u'\u03C9')]

    TOTAL_INDEX = {'na':0,'e':1,'E':2,'S':3,'h':4,'H':5,'t':6,'T':7,'g':8,'G':9,'B':10,'i':11,'I':12}
    STAT_INDEX = {}

    ss_field = 'r%i_ss' % iIndex
    aa_field = 'r%i_aa' % iIndex

    stat_attributes = ['r%i_%s'%(iIndex, f[0]) for f in stat_attributes_base]
    fieldNames      = ['r%i_%s'%(iIndex, f[0]) for f in stat_attributes_base]
    fieldNames.append(ss_field)

    for i in range(len(stat_attributes)):
        STAT_INDEX[stat_attributes[i]] = i

    # get search from session
    search = request.session['search']
    local_query_filter = search.querySet().filter

    peptides = {}

    total = {
                #total
                'total':0,
                # attributes with just sums
                'counts':[['na',0],['e',0],['E',0],['S',0],['h',0],['H',0],['t',0],['T',0],['g',0],['G',0],['B',0],['i',0],['I',0]],
                # attributes with stats
                'stats':[['L1',[]],['L2',[]],['L3',[]],['L4',[]],['L5',[]],['a1',[]],['a2',[]],['a3',[]],['a4',[]],['a5',[]],['a6',[]],['a7',[]],['ome',[]]]
            }

    #iterate through the aa_choices
    for code,long_code in AA_CHOICES:
        #create data structure
        peptide = {
                'longCode':long_code,
                #total
                'total':0,
                # attributes with just sums
                'counts':[['na',0],['e',0],['E',0],['S',0],['h',0],['H',0],['t',0],['T',0],['g',0],['G',0],['B',0],['i',0],['I',0]],
                # attributes with stats
                'stats':[['L1',[]],['L2',[]],['L3',[]],['L4',[]],['L5',[]],['a1',[]],['a2',[]],['a3',[]],['a4',[]],['a5',[]],['a6',[]],['a7',[]],['ome',[]]]
            }

        #query segments matching this AA with just the fields we want to perform calcuations on
        residueData = local_query_filter(**{aa_field:code}).values(*fieldNames)
        peptide['total'] = local_query_filter(**{aa_field:code}).count()
        total['total'] += peptide['total']

        #iterate through all the segment data
        for data in residueData:

            #calculate values
            if data[ss_field] in (' ','-'):
                peptide['counts'][TOTAL_INDEX['na']][1] += 1
                total['counts'][TOTAL_INDEX['na']][1] += 1
            else:
                peptide['counts'][TOTAL_INDEX[data[ss_field]]][1] += 1
                total['counts'][TOTAL_INDEX[data[ss_field]]][1] += 1

            #store all values for attributes into arrays
            for key in stat_attributes:
                peptide['stats'][STAT_INDEX[key]][1].append(data[key])
                total['stats'][STAT_INDEX[key]][1].append(data[key])

        #calculate statistics
        for attribute in stat_attributes:
            list = peptide['stats'][STAT_INDEX[attribute]][1]
            peptide['stats'][STAT_INDEX[attribute]][1] = calculate_stats(list, attribute)

        peptides[code] = peptide


    #calculate statistics for totals
    for attribute in stat_attributes:
        list = total['stats'][STAT_INDEX[attribute]][1]
        total['stats'][STAT_INDEX[attribute]][1] = calculate_stats(list, attribute)

    return render_to_response('stats.html', {
        'attributes': stat_attributes_base,
        'peptides':peptides,
        'total':total
    }, context_instance=RequestContext(request))


def calculate_stats(_list, attribute):
    """
    Calculates stats for a list of residues and saves it in the dictionary
    passed in.  The expected structure is defined  in the above view handler
    for statistics.
    """

    local_pow = pow

    list_len = len(_list)
    if list_len > 1:
        if attribute in angles:
            mean, stdev = getCircularStats(_list, list_len)
            range_min = 0
            range_max = 0
        else:
            #Average/mean
            mean = sum(_list)/list_len

            #now that we have mean calculate standard deviation
            stdev = math.sqrt(
                sum([
                    local_pow(value - mean, 2)
                    for value in _list
                ])/(list_len - 1)
            )

            range_min = '%+.3f' % (min(_list) - mean)
            range_max = '%+.3f' % (max(_list) - mean)

    # if theres only 1 item then the stats are simpler to calculate
    elif list_len == 1:
        mean = _list[0]
        std_dev = 0
        range_min = 0
        range_max = 0

    else:
        mean = 0
        stdev = 0
        range_min = 0
        range_max = 0

    return {'mean':mean,'std':stdev,'min':range_min, 'max':range_max}