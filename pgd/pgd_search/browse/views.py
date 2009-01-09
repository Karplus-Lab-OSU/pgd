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
    rStart = 50;
    gStart = 180;
    bStart = 50;
    rInterval = round((255-rStart)/stop)
    gInterval = round((255-gStart)/stop)
    bInterval = round((200-bStart)/stop)
    for i in range(stop):
        red = '%s' % hex(int(i*rInterval+rStart))[2:] if i*rInterval+rStart > 9 else '%02d' % int(hex(int(i*rInterval+rStart))[2:])
        green = '%s' % hex(int(i*gInterval+gStart))[2:] if i*gInterval+gStart > 9 else '%02d' % int(hex(int(i*gInterval+gStart))[2:])
        blue = '%s' % hex(int(i*bInterval+bStart))[2:] if i*bInterval+bStart > 9 else '%02d' % int(hex(int(i*bInterval+bStart))[2:])
        colors['.index%i, .index-%i' % (i,i)] = (red,green,blue)

    return render_to_response('browse.html', {
        'segments': paginatedSegments,
        'segmentSlice':'%i:%i'%(iIndex+start,iIndex+stop),
        'iIndex':iIndex,
        'pageStart':(page-1)*5,
        'indexColors':colors        
    }, context_instance=RequestContext(request))
