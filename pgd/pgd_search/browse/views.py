from django.shortcuts import render_to_response
from django.template import RequestContext
from django.core.paginator import Paginator, InvalidPage, EmptyPage
import math

from pgd_search.models import searchSettings

"""
display search results in tabular form
"""
def browse(request):

    # get search from session
    search = request.session['search']
    segments = search.querySet()
    paginator = Paginator(segments, 5) # Show 10 segments per page

    # Make sure page request is an int. If not, deliver first page.
    try:
        page = int(request.GET.get('page', '1'))
    except ValueError:
        page = 1

    # If page request (9999) is out of range, deliver last page of results.
    try:
        paginatedSegments = paginator.page(page)
    except (EmptyPage, InvalidPage):
        paginatedSegments = paginator.page(paginator.num_pages)

    #generate iValues
    iIndex = int(math.ceil(searchSettings.segmentSize/2.0)-1)
    start = 0 - (search.residueCount-1) / 2
    stop  = int(math.ceil((search.residueCount-1) / 2.0))+1
    iValues = [
            # Generates a series of tuples (<value>,<signed string of value>); zero is 'i'
            '%+i'%i if i else 'i' for i in range(start,stop)
    ]

    return render_to_response('browse.html', {
        'segments': paginatedSegments,
        'iValues' : iValues,
        'segmentSlice':'%i:%i'%(iIndex+start,iIndex+stop),
        'pageStart':(page-1)*5
    }, context_instance=RequestContext(request))
