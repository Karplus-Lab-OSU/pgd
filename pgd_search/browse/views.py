from django.conf import settings
from django.shortcuts import render_to_response
from django.template import RequestContext
from django.core.paginator import Paginator, InvalidPage, EmptyPage
import math
import pickle
from pgd_search.views import settings_processor

"""
display search results in tabular form
"""
def browse(request):
    #saving globals as locals for speed
    lhex = hex
    lint = int

    # get search from session
    search = pickle.loads(request.session['search'])

    #generate iValues
    start = 0 - (search.segmentLength-1) / 2
    stop  = lint(math.ceil((search.segmentLength-1) / 2.0))+1
    iIndex = lint(math.ceil(settings.SEGMENT_SIZE/2.0)-1)

    #paginate
    segments = search.querySet().order_by('protein')
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
        page = paginator.num_pages
        paginatedSegments = paginator.page(page)

    #generate a list of pages to display in the pagination bar
    pages = get_page_list(paginator, page)


    # use ranges for RGB to introduce colorful steps
    colors = [] 
    colorstop=settings.SEGMENT_SIZE/2
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
        'indexes':range(start,stop),
        'start':start,
        'stop':stop-1,
        'iIndex':iIndex,
        'pageStart':(page-1)*25,
        'indexColors':colors,
        'pages':pages,
        'segmentLength':search.segmentLength

    }, context_instance=RequestContext(request, processors=[settings_processor]))


def get_page_list(paginator, page=1):
    """
    Generate a list of pages used choose from to jump quickly to a page

    This generates a list that shows:
        * if near the start/end, up to 10 pages from the start/end
        * if in the middle first two and last two pages in addition to the
          +/- 5 from the current page
       * if num_pages<10 - only one list of pages is shown
       * if num_pages==11 then the list is statically generated because this
         list size results in wierd results for the standard logic that generates
         the lists.
    """

    if paginator.num_pages < 11:
        # special case: only one list of values
        pages = (range(1,paginator.num_pages+1),)
    elif paginator.num_pages == 11:
        # special case: lists are static
        pages = ([1,2,3,4,5,6,7,8,9,10], None, [11])
    else:
        # normal case
        start = [i for i in range(1, 11 if page < 8 and page < paginator.num_pages-6 else 3)]
        middle = [i for i in range(page-5,page+5)] if page > 7 and page < paginator.num_pages-6 else None
        end = [i for i in range(paginator.num_pages-(1 if page < paginator.num_pages-6 else 9), paginator.num_pages+1)]
        pages = (start, middle, end)

    return pages