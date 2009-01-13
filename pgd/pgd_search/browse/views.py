from django.shortcuts import render_to_response
from django.template import RequestContext
from django.core.paginator import Paginator, InvalidPage, EmptyPage
import math

from pgd_search.models import searchSettings

"""
display search results in tabular form
"""
def browse(request):

    #saving globals as locals for speed
    lhex = hex
    lint = int

    # get search from session
    search = request.session['search']
    segments = search.querySet()
    paginator = Paginator(segments, 25) # Show 25 segments per page

    # Make sure page request is an int. If not, deliver first page.
    try:
        page = lint(request.GET.get('page', '1'))
    except ValueError:
        page = 1

    # If page request (9999) is out of range, deliver last page of results.
    try:
        paginatedSegments = paginator.page(page)
    except (EmptyPage, InvalidPage):
        paginatedSegments = paginator.page(paginator.num_pages)

    #generate iValues
    iIndex = lint(math.ceil(searchSettings.segmentSize/2.0)-1)
    start = 0 - (search.segmentLength-1) / 2
    stop  = lint(math.ceil((search.segmentLength-1) / 2.0))+1

    # use ranges for RGB to introduce colorful steps
    colors = [] 
    colorstop=searchSettings.segmentSize/2
    rRange = (75,245);
    gRange = (180,245);
    bRange = (240,245);
    rInterval = round((rRange[1]-rRange[0])/colorstop)
    gInterval = round((gRange[1]-gRange[0])/colorstop)
    bInterval = round((bRange[1]-bRange[0])/colorstop)

    for i in range(colorstop-stop+1, colorstop+1):
        red = '%s' % lhex(lint(i*rInterval+rRange[0]))[2:] if i*rInterval+rRange[0] > 9 else '%02d' % lint(lhex(lint(i*rInterval+rRange[0]))[2:])
        green = '%s' % lhex(lint(i*gInterval+gRange[0]))[2:] if i*gInterval+gRange[0] > 9 else '%02d' % lint(lhex(lint(i*gInterval+gRange[0]))[2:])
        blue = '%s' % lhex(lint(i*bInterval+bRange[0]))[2:] if i*bInterval+bRange[0] > 9 else '%02d' % lint(lhex(lint(i*bInterval+bRange[0]))[2:])
        colors.append((red,green,blue))

    return render_to_response('browse.html', {
        'segments': paginatedSegments,
        'segmentSlice':'%i:%i'%(iIndex+start,iIndex+stop),
        'iIndex':iIndex,
        'pageStart':(page-1)*5,
        'indexColors':colors
    }, context_instance=RequestContext(request))
