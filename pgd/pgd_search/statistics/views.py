from django.shortcuts import render_to_response
from django.template import RequestContext
from django.conf import settings
import math

from constants import AA_CHOICES
from pgd_search.models import Search, Segment, iIndex
from pgd_search.plot.ConfDistFuncs import getCircularStats

"""
display statistics about the search
"""
def searchStatistics(request):

    #store globals locally for speed optimization
    local_len = len
    local_sum = sum
    local_pow = pow
    local_min = min
    local_max = max

    ANGLES_BASE = ('ome', 'phi', 'psi', 'chi', 'zeta')
    angles = ['r%i_%s' %(iIndex, angle) for angle in ANGLES_BASE]

    stat_attributes_base = ['L1','L2','L3','L4','L5','a1','a2','a3','a4','a5','a6','a7','ome']
    TOTAL_INDEX = {'na':0,'e':1,'E':2,'S':3,'h':4,'H':5,'t':6,'T':7,'g':8,'G':9,'B':10,'i':11,'I':12}
    STAT_INDEX = {}

    ss_field = 'r%i_ss' % iIndex
    aa_field = 'r%i_aa' % iIndex

    stat_attributes = ['r%i_%s'%(iIndex, f) for f in stat_attributes_base]
    fieldNames      = ['r%i_%s'%(iIndex, f) for f in stat_attributes_base]
    fieldNames.append(ss_field)

    for i in range(local_len(stat_attributes)):
        STAT_INDEX[stat_attributes[i]] = i

    # get search from session
    search = request.session['search']
    local_query_filter = search.querySet().filter

    peptides = {}
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

        #iterate through all the segment data
        for data in residueData:

            #calculate values
            if data[ss_field] == ' ':
                peptide['counts'][TOTAL_INDEX['na']][1] += 1
            else:
                peptide['counts'][TOTAL_INDEX[data[ss_field]]][1] += 1

            #store all values for attributes into arrays
            for key in stat_attributes:
                peptide['stats'][STAT_INDEX[key]][1].append(data[key])

        #calculate statistics
        for attribute in stat_attributes:
            list = peptide['stats'][STAT_INDEX[attribute]][1]

            list_len = local_len(list)
            if list_len > 1:
                if attribute in angles:
                    mean, std_dev = getCircularStats(list, local_len(list))
                    range_min = 0
                    range_max = 0
                else:
                    #Average/mean
                    mean = local_sum(list)/list_len

                    #now that we have mean calculate standard deviation
                    stdev = math.sqrt(
                        local_sum([
                            local_pow(value - mean, 2)
                            for value in list
                        ])/(list_len - 1)
                    )

                    range_min = '%+.3f' % (local_min(list) - mean)
                    range_max = '%+.3f' % (local_max(list) - mean)

            # if theres only 1 item then the stats are simpler to calculate
            elif list_len == 1:
                mean = list[0]
                std_dev = 0
                range_min = 0
                range_max = 0

            else:
                mean = 0
                stdev = 0
                range_min = 0
                range_max = 0

            peptide['stats'][STAT_INDEX[attribute]][1] = {'mean':mean,'std':stdev,'min':range_min, 'max':range_max}

        peptides[code] = peptide

    return render_to_response('stats.html', {
        'attributes': stat_attributes_base,
        'peptides':peptides
    }, context_instance=RequestContext(request))