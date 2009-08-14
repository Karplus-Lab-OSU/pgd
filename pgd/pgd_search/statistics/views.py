from django.shortcuts import render_to_response
from django.template import RequestContext
from django.conf import settings
import math

from pgd_constants import AA_CHOICES
from pgd_search.models import Search, Segment, iIndex
from pgd_search.plot.ConfDistFuncs import getCircularStats
from pgd_search.views import settings_processor

stat_attributes_base = [('L1',u'C<sup>-1</sup>N'),
                        ('L2',u'NC<sup>&alpha;</sup>'),
                        ('L3',u'C<sup>&alpha;</sup>C<sup>&beta;</sup>'),
                        ('L4',u'C<sup>&alpha;</sup>C'),
                        ('L5','CO'),
                        ('a1',u'C<sup>-1</sup>NC<sup>&alpha;</sup>'),
                        ('a2',u'NC<sup>&alpha;</sup>C<sup>&beta;</sup>'),
                        ('a3',u'NC<sup>&alpha;</sup>C'),
                        ('a4',u'C<sup>&beta;</sup>C<sup>&alpha;</sup>C'),
                        ('a5',u'C<sup>&alpha;</sup>CO'),
                        ('a6',u'C<sup>&alpha;</sup>CN<sup>+1</sup>'),
                        ('a7',u'OCN<sup>+1</sup>'),
                        ('ome',u'&omega;')]

ANGLES_BASE = ('ome', 'phi', 'psi', 'chi', 'zeta')
angles = ['r%i_%s' %(iIndex, angle) for angle in ANGLES_BASE]


"""
display statistics about the search
"""
def searchStatistics(request):

    # get search from session
    search = request.session['search']

    peptides, total = calculate_statistics(search.querySet())

    return render_to_response('stats.html', {
        'attributes': stat_attributes_base,
        'peptides':peptides,
        'total':total
    }, context_instance=RequestContext(request, processors=[settings_processor]))



def calculate_statistics(queryset):

    TOTAL_INDEX = {'na':0,'E':1,'S':2,'H':3,'T':4,'G':5,'B':6,'I':7}
    STAT_INDEX = {}

    ss_field = 'r%i_ss' % iIndex
    aa_field = 'r%i_aa' % iIndex

    stat_attributes = ['r%i_%s'%(iIndex, f[0]) for f in stat_attributes_base]
    fieldNames      = ['r%i_%s'%(iIndex, f[0]) for f in stat_attributes_base]
    fieldNames.append(ss_field)

    for i in range(len(stat_attributes)):
        STAT_INDEX[stat_attributes[i]] = i


    local_query_filter = queryset.filter

    peptides = {}

    total = {
                #total
                'total':0,
                # attributes with just sums
                'counts':[['na',0],['E',0],['S',0],['H',0],['T',0],['G',0],['B',0],['I',0]],
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
                'counts':[['na',0],['E',0],['S',0],['H',0],['T',0],['G',0],['B',0],['I',0]],
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
                if not data[key] in (999.9, 0):
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

    return peptides, total


def calculate_stats(_list, attribute):
    """
    Calculates stats for a list of residues and saves it in the dictionary
    passed in.  The expected structure is defined  in the above view handler
    for statistics.
    """

    #store globals locally for speed optimization
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
        stdev = 0
        range_min = 0
        range_max = 0

    else:
        mean = 0
        stdev = 0
        range_min = 0
        range_max = 0

    return {'mean':mean,'std':stdev,'min':range_min, 'max':range_max}
