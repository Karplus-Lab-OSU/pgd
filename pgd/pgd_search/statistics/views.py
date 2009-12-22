import math
import simplejson
import time
from threading import Thread

from django.conf import settings
from django.db.models import Avg, Max, Min, Count, StdDev
from django.http import HttpResponse
from django.shortcuts import render_to_response
from django.template import RequestContext

from pgd_constants import AA_CHOICES, SS_CHOICES
from pgd_search.models import Search, Segment
from pgd_search.statistics.aggregates import *
from pgd_search.statistics.directional_stddev import *
from pgd_search.views import settings_processor



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
                        ('ome',u'&omega;')]

FIELDS_BASE = ('L1','L2','L3','L4','L5','a1','a2','a3','a4','a5','a6','a7')
ANGLES_BASE = ('ome', ) #,'phi')#, 'psi', 'chi', 'zeta')


def search_statistics(request):
    """
    display statistics about the search
    """
    return render_to_response('stats.html', {
        'aa_types': AA_CHOICES,
        'ss_types': SS_CHOICES,
        'fields':FIELDS_BASE,
        'angles':ANGLES_BASE,
        'stat_attributes':stat_attributes
    }, context_instance=RequestContext(request, processors=[settings_processor]))


def search_statistics_data(request):
    """
    returns ajax'ified statistics data for the current search
    """
    search = request.session['search']
    try:
        stats = calculate_statistics_a(search.querySet())
    except Exception, e:
        print 'exception', e
        import traceback, sys
        exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        print "*** print_tb:"
        traceback.print_tb(exceptionTraceback, limit=10, file=sys.stdout)

        raise e
    return HttpResponse(simplejson.dumps(stats))


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

    # statistics for fields - build a list of aggregations to annote
    # against all fields that calculate avg/std normally
    field_annotations = {}
    for f in FIELDS_BASE:
        field_annotations['avg_%s' % f] = Avg(field_prefix % f)
        field_annotations['stddev_%s' % f] = StdDev(field_prefix % f)
        field_annotations['min_%s' % f] = Min(field_prefix % f)
        field_annotations['max_%s' % f] = Max(field_prefix % f)
    field_stats = list(queryset.values(aa_field).annotate(**field_annotations) \
                                                .order_by(aa_field))
    field_stats_totals = queryset.aggregate(**field_annotations)

    now = time.time()
    print '  -fields: %s (%s)' % (now-start, now-last)
    last = now

    # statistics for angles - must be done in two passes because there is no
    # apparent way to do standard devation in a single pass
    angle_annotations = {}
    for f in ANGLES_BASE:
        angle_annotations['avg_%s' % f] = DirectionalAvg(field_prefix % f)
        angle_annotations['min_%s' % f] = Min(field_prefix % f)
        angle_annotations['max_%s' % f] = Max(field_prefix % f)
    angles_stats = list(queryset.values(aa_field).annotate(**angle_annotations) \
                                                    .order_by(aa_field))
    angles_stats_totals = queryset.aggregate(**angle_annotations)
    now = time.time()
    print '  -angles: %s (%s)' % (now-start, now-last)
    last = now

    angles_stddev = []
    for res in angles_stats:
        aa = res[aa_field]
        angle_annotations = {}
        for f in ANGLES_BASE:
            field = field_prefix % f
            avg = res['avg_%s' % f]
            angle_annotations['stddev_%s' % f] = DirectionalStdDev(field, avg=avg)
        q = queryset.filter(**{aa_field:aa}).aggregate(**angle_annotations)
        q[field] = aa
        angles_stddev.append(q)
    now = time.time()
    print '  -angles stdev: %s (%s)' % (now-start, now-last)
    last = now

    angle_annotations = {}
    for f in ANGLES_BASE:
        field = field_prefix % f
        avg = angles_stats_totals['avg_%s' % f]
        angle_annotations['stddev_%s' % f] = DirectionalStdDev(field, avg=avg)
    angles_stddev_totals = queryset.aggregate(**angle_annotations)

    now = time.time()
    print '  -angles total stdev: %s (%s)' % (now-start, now-last)
    last = now


    # ss/aa counts and totals
    ss_counts = list(queryset.values(aa_field, ss_field).annotate(ss_count=Count(ss_field)).order_by(aa_field))
    ss_totals = list(queryset.values(ss_field).annotate(ss_count=Count(ss_field)))
    aa_totals = list(queryset.values(aa_field).annotate(aa_count=Count(aa_field)).order_by(aa_field))
    now = time.time()
    print '  -counts: %s (%s)' % (now-start, now-last)
    last = now

    stats = {
        'index':4,
        'fields': field_stats,                       # list of dictionaries, list by aa
        'fields_totals':field_stats_totals,          # dict of agg per field
        'angles':angles_stats,                       # list of dictionaries, list by aa
        'angles_stddev':angles_stddev,               # list of dictionaries, list by aa
        'angles_totals':angles_stats_totals,         # dictionary of agg per field
        'angles_stddev_totals':angles_stddev_totals, # dictionary of agg per field
        'aa_totals':aa_totals,
        'ss_counts':ss_counts,                       # list of dictionaries, list by AA/SS
        'ss_totals':ss_totals,                        # list of dictionaries,  list by SS
        'total':queryset.count()
    }

    now = time.time()
    print '  -stats: %s (%s)' % (now-start, now-last)
    last = now
    end = time.time()
    print 'Search Statistics Data in seconds: ', end-start

    return stats


def calculate_statistics_a(queryset, iIndex=0):
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

    # get all stats fields - this incldues both regular fields plus dihedral angles.  This does not include
    # field totals yet but it will some time in the future.  Start this rthread immediately because this 
    # query is 3x as long as the other queries
    dsq_thread = ListQueryThread(DirectionalStatisticsQuery(ANGLES_BASE, FIELDS_BASE, field_prefix, queryset))
    dsq_thread.start()

    # ss/aa counts and totals
    ss_counts_thread = ListQueryThread(queryset.values(aa_field, ss_field).annotate(ss_count=Count(ss_field)))
    ss_totals_thread = ListQueryThread(queryset.values(ss_field).annotate(ss_count=Count(ss_field)))
    aa_totals_thread = ListQueryThread(queryset.values(aa_field).annotate(aa_count=Count(aa_field)))

    # start all remtaining threads
    ss_counts_thread.start()
    ss_totals_thread.start()
    aa_totals_thread.start()

    # wait for each thread independantly by joining with the thread. 
    # as threads finish they will return.  If a thread is already done
    # it will return.  When all threads have returned they are all done 
    # and its safe to ask them for the results
    dsq_thread.join()
    ss_counts_thread.join()
    ss_totals_thread.join()
    aa_totals_thread.join()

    angles_stats = dsq_thread.results
    ss_counts = ss_counts_thread.results
    ss_totals = ss_totals_thread.results
    aa_totals = aa_totals_thread.results

    stats = {
        'index':4,
        'fields': angles_stats,                       # list of dictionaries, list by aa
        'angles_totals':[],         # dictionary of agg per field
        'angles_stddev_totals':[], # dictionary of agg per field
        'aa_totals':aa_totals,
        'ss_counts':ss_counts,                       # list of dictionaries, list by AA/SS
        'ss_totals':ss_totals,                        # list of dictionaries,  list by SS
        'total':queryset.count()
    }

    now = time.time()
    print '  -stats: %s (%s)' % (now-start, now-last)
    last = now
    end = time.time()
    print 'Search Statistics Data in seconds: ', end-start

    return stats


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
        print '  -counts: f=%s, t=%s' % (end, end-start)
