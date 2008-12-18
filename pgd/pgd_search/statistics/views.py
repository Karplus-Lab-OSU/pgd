from django.shortcuts import render_to_response
from django.conf import settings
from statlib import stats

from constants import AA_CHOICES
from pgd_search.models import Search, Segment


"""
display statistics about the search
"""
def searchStatistics(request):

    stat_attributes = ['L1','L2','L3','L4','L5','a1','a2','a3','a4','a5','a6','a7']
    TOTAL_INDEX = {'na':0,'e':1,'E':2,'S':3,'h':4,'H':5,'t':6,'T':7,'g':8,'G':9,'B':10,'i':11,'I':12}
    STAT_INDEX = {'L1':0,'L2':1,'L3':2,'L4':3,'L5':4,'a1':5,'a2':6,'a3':7,'a4':8,'a5':9,'a6':10,'a7':11}

    # get search from session
    search = Search()

    # parse search into djangoQuery
    #searchQuery = queryParser(search)
    searchQuery = Segment.objects.all()[:500]#temp replacement for testing

    peptides = {}
    for code,long_code in AA_CHOICES:
        peptides[code] = {
                'longCode':long_code,
                #total
                'total':0,
                # attributes with just sums
                'counts':[['na',0],['e',0],['E',0],['S',0],['h',0],['H',0],['t',0],['T',0],['g',0],['G',0],['B',0],['i',0],['I',0]],
                # attributes with stats
                'stats':[['L1',[]],['L2',[]],['L3',[]],['L4',[]],['L5',[]],['a1',[]],['a2',[]],['a3',[]],['a4',[]],['a5',[]],['a6',[]],['a7',[]]]
            }

    #iterate through all the segments
    for segment in searchQuery:

        #get the i protein
        segment = segment.residues[4]

        #get the peptide
        peptide = peptides[segment.aa]

        #calculate values
        peptide['total'] += 1
        if segment.ss == ' ':
            peptide['counts'][TOTAL_INDEX['na']][1] += 1
        else:
            peptide['counts'][TOTAL_INDEX[segment.ss]][1] += 1

        #store all values for attributes into arrays
        for key in stat_attributes:
            peptide['stats'][STAT_INDEX[key]][1].append(segment.__dict__[key])


    # iterate through all the datastructures calculating statistics
    for key in peptides:
        peptide = peptides[key]
        for attribute in stat_attributes:
            list = peptide['stats'][STAT_INDEX[attribute]][1]
            list_len = len(list)
            if list_len > 1:
                mean = stats.mean(list)
                #now that we have mean calculate standard deviation
                stdev = stats.stdev(list)
                range_min = '%+.3f' % (min(list) - mean)
                range_max = '%+.3f' % (max(list) - mean)

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


    return render_to_response('stats.html', {
        'SITE_ROOT': settings.SITE_ROOT,
        'MEDIA_URL': settings.MEDIA_URL,
        'attributes': stat_attributes,
        'peptides':peptides
    })