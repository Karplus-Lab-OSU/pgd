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
    start = 0 - (search.segmentLength-1) / 2
    stop  = int(math.ceil((search.segmentLength-1) / 2.0))+1

    # use ranges for RGB to introduce colorful steps
    colors = {} 
    rRange = (75,200);
    gRange = (180,240);
    bRange = (240,255);
    rInterval = round((rRange[1]-rRange[0])/stop)
    gInterval = round((gRange[1]-gRange[0])/stop)
    bInterval = round((bRange[1]-bRange[0])/stop)
    for i in range(stop):

        red = '%s' % hex(int(i*rInterval+rRange[0]))[2:] if i*rInterval+rRange[0] > 9 else '%02d' % int(hex(int(i*rInterval+rRange[0]))[2:])
        green = '%s' % hex(int(i*gInterval+gRange[0]))[2:] if i*gInterval+gRange[0] > 9 else '%02d' % int(hex(int(i*gInterval+gRange[0]))[2:])
        blue = '%s' % hex(int(i*bInterval+bRange[0]))[2:] if i*bInterval+bRange[0] > 9 else '%02d' % int(hex(int(i*bInterval+bRange[0]))[2:])
        colors['.index%i, .index-%i' % (i,i)] = (red,green,blue)

    return render_to_response('browse.html', {
        'segments': paginatedSegments,
        'segmentSlice':'%i:%i'%(iIndex+start,iIndex+stop),
        'iIndex':iIndex,
        'pageStart':(page-1)*5,
        'indexColors':colors        
    }, context_instance=RequestContext(request))
