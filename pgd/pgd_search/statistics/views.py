import math
import time
import simplejson

from django.conf import settings
from django.db.models import Avg, Max, Min, Count, StdDev
from django.http import HttpResponse
from django.shortcuts import render_to_response
from django.template import RequestContext

from pgd_constants import AA_CHOICES, SS_CHOICES
from pgd_search.models import Search, Segment, iIndex
from pgd_search.statistics.aggregates import *
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

FIELDS_BASE = ('L1','L2','L3','L4','L5','a1','a2','a3','a4','a5','a6','a7')
ANGLES_BASE = ('ome', ) #,'phi')#, 'psi', 'chi', 'zeta')
angles = ['r%i_%s' %(iIndex, angle) for angle in ANGLES_BASE]


def search_statistics(request):
    """
    display statistics about the search
    """
    return render_to_response('stats.html', {
        'aa_types': AA_CHOICES,
        'ss_types': SS_CHOICES,
        'fields':FIELDS_BASE,
        'angles':ANGLES_BASE
    }, context_instance=RequestContext(request, processors=[settings_processor]))


def search_statistics_data(request):
    """
    returns ajax'ified statistics data for the current search
    """
    search = request.session['search']
    stats = calculate_statistics(search.querySet())

    return HttpResponse(simplejson.dumps(stats))


def calculate_statistics(queryset):
    """
    Calculates statistics across most fields stored in the database.  This
    function uses SQL aggregate functions to perform calculations in place
    on the database.  This requires multiple queries but ultimately performs
    and scales better with large datasets.
    """

    start = time.time()
    ss_field = 'r%i_ss' % iIndex
    aa_field = 'r%i_aa' % iIndex

    # statistics for fields - build a list of aggregations to annote
    # against all fields that calculate avg/std normally
    field_annotations = {}
    for f in FIELDS_BASE:
        field_annotations['avg_%s' % f] = Avg('r4_%s' % f)
        field_annotations['stddev_%s' % f] = StdDev('r4_%s' % f)
        field_annotations['min_%s' % f] = Min('r4_%s' % f)
        field_annotations['max_%s' % f] = Max('r4_%s' % f)
    field_stats = queryset.values(aa_field).annotate(**field_annotations) \
                                                .order_by(aa_field)
    field_stats_totals = queryset.aggregate(**field_annotations)

    # statistics for angles - must be done in two passes because there is no
    # apparent way to do standard devation in a single pass
    angle_annotations = {}
    for f in ANGLES_BASE:
        angle_annotations['avg_%s' % f] = DirectionalAvg('r4_%s' % f)
        angle_annotations['min_%s' % f] = Min('r4_%s' % f)
        angle_annotations['max_%s' % f] = Max('r4_%s' % f)
    angles_stats = queryset.values(aa_field).annotate(**angle_annotations) \
                                                    .order_by(aa_field)

    angles_stats_totals = queryset.aggregate(**angle_annotations)

    angles_stddev = []
    for res in angles_stats:
        aa = res[aa_field]
        angle_annotations = {}
        for f in ANGLES_BASE:
            field = 'r4_%s' % f
            avg = res['avg_%s' % f]
            angle_annotations['stddev_%s' % f] = DirectionalStdDev(field, avg=avg)
        angles_stddev.append(queryset.filter(**{aa_field:aa}).values(aa_field).annotate(**angle_annotations))

    angle_annotations = {}
    for f in ANGLES_BASE:
        field = 'r4_%s' % f
        avg = angles_stats_totals['avg_%s' % f]
        angle_annotations['stddev_%s' % f] = DirectionalStdDev(field, avg=avg)
    angles_stddev_totals = queryset.aggregate(**angle_annotations)

    # ss/aa counts and totals
    ss_counts = queryset.values(aa_field, ss_field).annotate(ss_count=Count(ss_field)).order_by(aa_field)
    ss_totals = queryset.values(ss_field).annotate(ss_count=Count(ss_field))
    aa_totals = queryset.values(aa_field).annotate(aa_count=Count(aa_field)).order_by(aa_field)

    end = time.time()
    print 'Search Statistics Data in seconds: ', end-start

    return {
        'index':4,
        'fields': list(field_stats),                       # list of dictionaries, list by aa
        'fields_totals':field_stats_totals,          # dict of agg per field
        'angles':list(angles_stats),                       # list of dictionaries, list by aa
        'angles_stddev':[list(s) for s in angles_stddev],               # list of dictionaries, list by aa
        'angles_totals':angles_stats_totals,         # dictionary of agg per field
        'angles_stddev_totals':angles_stddev_totals, # dictionary of agg per field
        'aa_totals':list(aa_totals),
        'ss_counts':list(ss_counts),                       # list of dictionaries, list by AA/SS
        'ss_totals':list(ss_totals),                        # list of dictionaries,  list by SS
        'total':queryset.count()
    }