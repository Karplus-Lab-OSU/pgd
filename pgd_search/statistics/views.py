import json
import time
import pickle
from threading import Thread

from django.db.models import Count
from django.http import HttpResponse
from django.shortcuts import render_to_response
from django.template import RequestContext
from pgd_constants import AA_CHOICES, SS_CHOICES, AA_CHOICES_DICT
from pgd_search.statistics.aggregates import *
from pgd_search.statistics.directional_stddev import *
from pgd_search.statistics.form import StatsForm
from pgd_search.views import settings_processor
from pgd_splicer.sidechain import bond_angles_string_dict, bond_lengths_string_dict

stat_attributes = [('L1',u'C<sup>-1</sup>N'),
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
                        ('ome',u'&omega;'),
                        ("omep",u'&omega;P')]

FIELDS_BASE = ('L1','L2','L3','L4','L5','a1','a2','a3','a4','a5','a6','a7')
ANGLES_BASE = ('ome', 'omep') #,'phi')#, 'psi', 'chi', 'zeta')


LENGTHS = [('L1',u'C<sup>-1</sup>N'),
        ('L2',u'NC<sup>&alpha;</sup>'),
        ('L3',u'C<sup>&alpha;</sup>C<sup>&beta;</sup>'),
        ('L4',u'C<sup>&alpha;</sup>C'),
        ('L5','CO')]


ANGLES = [('a1',u'C<sup>-1</sup>NC<sup>&alpha;</sup>'),
        ('a2',u'NC<sup>&alpha;</sup>C<sup>&beta;</sup>'),
        ('a3',u'NC<sup>&alpha;</sup>C'),
        ('a4',u'C<sup>&beta;</sup>C<sup>&alpha;</sup>C'),
        ('a5',u'C<sup>&alpha;</sup>CO'),
        ('a6',u'C<sup>&alpha;</sup>CN<sup>+1</sup>'),
        ('a7',u'OCN<sup>+1</sup>'),
        ('ome',u'&omega;'),
        ("omep",u'&omega;P')]


# rebuild bond angle list so they work with single letter codes
BOND_ANGLES = {}
BOND_LENGTHS = {}
for k, v in AA_CHOICES_DICT.items():
    v = v.upper()
    if v in bond_lengths_string_dict:
        BOND_LENGTHS[k] = [s.replace('-','_') for s in bond_lengths_string_dict[v]]
        BOND_ANGLES[k] = [s.replace('-','_') for s in bond_angles_string_dict[v]]


def search_statistics(request):
    """
    display statistics about the search
    """
    fields = dict(
        lengths = LENGTHS,
        angles = ANGLES,
        sca = BOND_ANGLES,
        scl = BOND_LENGTHS
    )
    
    return render_to_response('stats.html', {
        'form': StatsForm(),
        'aa_types': AA_CHOICES,
        'ss_types': SS_CHOICES,
        'fields':FIELDS_BASE,
        'stat_fields':fields,
        'angles':ANGLES_BASE,
        'stat_attributes':stat_attributes
    }, context_instance=RequestContext(request, processors=[settings_processor]))


def search_statistics_data(request):
    """
    returns ajax'ified statistics data for the current search
    """
    search = pickle.loads(request.session['search'])
    try:        
        index = int(request.GET['i']) if request.GET.has_key('i') else 0
        stats = calculate_statistics(search.querySet(), index)
    except Exception, e:
        print 'exception', e
        import traceback, sys
        exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        print "*** print_tb:"
        traceback.print_tb(exceptionTraceback, limit=10, file=sys.stdout)

        raise e
    return HttpResponse(json.dumps(stats))


def search_statistics_aa_data(request, aa):
    """
    returns ajax'ified statistics data for a single aa in the current search
    """
    search = pickle.loads(request.session['search'])
    try:        
        index = int(request.GET['i']) if request.GET.has_key('i') else 0
        stats = calculate_aa_statistics(search.querySet(), aa, index)
    except Exception, e:
        print 'exception', e
        import traceback, sys
        exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        print "*** print_tb:"
        traceback.print_tb(exceptionTraceback, limit=10, file=sys.stdout)

        raise e
    return HttpResponse(json.dumps(stats))


def calculate_statistics(queryset, iIndex=0):
    """
    Calculates statistics across most fields stored in the database.  This
    function uses SQL aggregate functions to perform calculations in place
    on the database.  This requires multiple queries but ultimately performs
    and scales better with large datasets.
    """
    start = time.time()
    last = start

    # get field prefix for this residue
    if iIndex == 0:
        prefix = ''
    elif iIndex < 0:
        prefix = ''.join(['prev__' for i in range(iIndex, 0)])
    else:
        prefix = ''.join(['next__' for i in range(iIndex)])
    field_prefix = '%s%%s' % prefix  

    ss_field = '%sss' % prefix
    aa_field = '%saa' % prefix

    # get all stats fields - this includes both regular fields plus dihedral
    # angles.  Start this thread immediately because this  query is 3x as long
    # as the other queries
    dsq_thread = ListQueryThread(DirectionalStatisticsQuery(ANGLES_BASE, FIELDS_BASE, field_prefix, queryset))
    dsq_thread.start()
    dsqt_thread = ListQueryThread(DirectionalStatisticsTotalQuery(ANGLES_BASE, FIELDS_BASE, field_prefix, queryset))
    dsqt_thread.start()
   
    # ss/aa counts and totals
    ss_counts_thread = ListQueryThread(queryset.values(aa_field, ss_field).annotate(ss_count=Count(ss_field)))
    ss_totals_thread = ListQueryThread(queryset.values(ss_field).annotate(ss_count=Count(ss_field)))
    aa_totals_thread = ListQueryThread(queryset.values(aa_field).order_by(aa_field).annotate(aa_count=Count(aa_field)))

    # start all remtaining threads
    ss_counts_thread.start()
    ss_totals_thread.start()
    aa_totals_thread.start()

    # wait for each thread independantly by joining with the thread. 
    # as threads finish they will return.  If a thread is already done
    # it will return.  When all threads have returned they are all done 
    # and its safe to ask them for the results
    dsq_thread.join()
    dsqt_thread.join()
    ss_counts_thread.join()
    ss_totals_thread.join()
    aa_totals_thread.join()

    field_stats = dsq_thread.results+dsqt_thread.results
    ss_counts = ss_counts_thread.results
    ss_totals = ss_totals_thread.results
    aa_totals = aa_totals_thread.results

    stats = {
        'prefix':prefix,
        'index':iIndex,
        'fields': field_stats,                       # list of dictionaries, list by aa
        'aa_totals':aa_totals,
        'ss_counts':ss_counts,                       # list of dictionaries, list by AA/SS
        'ss_totals':ss_totals,                        # list of dictionaries,  list by SS
        'total':queryset.count()
    }

    now = time.time()
    #print '  -stats: %s (%s)' % (now-start, now-last)
    last = now
    end = time.time()
    print 'Search Statistics Data in seconds: ', end-start

    return stats


def calculate_aa_statistics(queryset, aa, iIndex=0):
    """ Calculates statistics for a single AA type """
    # get field prefix for this residue
    if iIndex == 0:
        prefix = ''
    elif iIndex < 0:
        prefix = ''.join(['prev__' for i in range(iIndex, 0)])
    else:
        prefix = ''.join(['next__' for i in range(iIndex)])
    field_prefix = '%s%%s' % prefix  
    
    queryset = queryset.filter(**{'%saa'%prefix:aa})
    full_aa = AA_CHOICES_DICT[aa].upper()
    
    angles = []
    if aa in BOND_LENGTHS:
        fields = [str('sidechain_%s__%s' % (full_aa,f)) for f in BOND_ANGLES[aa]] + \
             [str('sidechain_%s__%s' % (full_aa,f)) for f in BOND_LENGTHS[aa]]
    else:
        fields = []
    
    dsq_thread = ListQueryThread(DirectionalStatisticsTotalQuery(angles, fields, field_prefix, queryset))
    dsq_thread.start()
    dsq_thread.join()
    results = dsq_thread.results[0]
    # XXX the Total query sets aa to total, but really its just a single AA
    results['aa'] = aa
    return results


class ListQueryThread(Thread):
    """
    Threaded class for running queries that return a list of the values
    This class is used to spin off queries into separate threads.  This
    allows stats view to be constrained only by the longest query.
    """
    def __init__(self, query):
        self.query = query
        Thread.__init__(self)

    def run(self):
        start = time.time()
        self.results = list(self.query)
        end = time.time()
        #print '  -counts: f=%s, t=%s' % (end, end-start)
