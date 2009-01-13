from django.http import HttpResponse
from django.template import RequestContext
from django.conf import settings
from django.shortcuts import render_to_response

from PlotForm import PlotForm
from ConfDistFuncs import *
from svg import *



"""
Renders a conformational distribution graph
@return: retusns an SVG instance.
"""
def drawGraph(request, xStart=-180.0, yStart=-180.0, xEnd=180.0, yEnd=180.0, attribute='Observations', xProperty='phi', yProperty='psi', reference=None, residue=None, xBin=10, yBin=10):
    
    #store functions locally for speed optimization
    local_len = len
    local_int = int

    svg = SVG()

    x = 65;
    y = 45;
    height = 360;
    width = 360;
    hashsize = 10

    #background
    svg.rect(x, y, width, height, 1, '#00fffff', '#222222');
    #svg.rect(0, 0, width+90, height+90, 1, '#00fffff');
    #svg.rect(x, y, width, height, 1, '#666666');

    #border
    svg.rect(x, y, width, height, 1, '#000000');

    #axis
    if xStart < 0 and xEnd > 0:
        xZero = (width/(xEnd-xStart)) * abs (xStart)
        svg.line( x+xZero, y, x+xZero, y+height, 1, '#666666');

    if yStart < 0 and xEnd > 0:
        yZero = height+y - (height/(yEnd-yStart)) * abs (yStart)
        svg.line( x, yZero, x+width, yZero, 1, '#666666');

    #hashes
    for i in range(9):
        hashx = x+(width/8.0)*i
        hashy = y+(height/8.0)*i
        svg.line( hashx, y+height, hashx, y+height+hashsize, 1, '#000000');
        svg.line( x, hashy, x-hashsize, hashy, 1, '#000000');

    #labels
    xstep = (xEnd - xStart) / 4
    ystep = (yEnd - yStart) / 4
    for i in range(5):
        xtext = xStart + xstep*i
        xtext = '%s' % local_int(xtext) if not xtext%1 else xtext
        xhash = x+(width/4)*i-(2.6*local_len(xtext))
        svg.text(xhash, y+height+hashsize*2+3, xtext,12)

        ytext = yEnd - ystep*i
        ytext = '%s' % local_int(ytext) if not ytext%1 else ytext
        yhash = y+(height/4)*i+4
        svg.text(x-15-(6*local_len(ytext)), yhash, ytext,12)

    #title text
    len1 = 220 - local_len(xProperty)*7/2 - local_len(xProperty)*7/2
    len2 = 182 - local_len(attribute)*7/2
    svg.text(len1,15, 'Plot of %s vs. %s' % (xProperty,yProperty), 12)
    svg.text(len2,35, 'Shading Based Off of %s' % attribute, 12)

    cdp = ConfDistPlot(
            360,            #height
            360,            #width
            0,              #Xpadding
            0,              #Ypadding
            x,              #Xoffset
            y,              #Yoffset
            xStart,         #Xstart
            xEnd,           #Xend
            yStart,         #Ystart
            yEnd,           #Yend
            xBin,           #Xbin
            yBin,           #Ybin
            xProperty,      #X property
            yProperty,      #Y property
            attribute,      #property
            residue,         #residue Index
            #reference
            request.session['search'].querySet()
    )

    boxes = cdp.Plot()
    return (svg,boxes)

def RGBTuple(rgbString):
    sub = rgbString[-6:]
    red = int(sub[:2],16)/255.0
    green = int(sub[2:4],16)/255.0
    blue = int(sub[4:], 16)/255.0
    return (red,green,blue)

def line(input, context):
    context.move_to(input.x+.5, input.y+.5)
    context.line_to(input.x1+.5, input.y1+.5)
    context.set_line_width(input.stroke)
    r,g,b = RGBTuple(input.color)
    context.set_source_rgba(r,g,b,1)
    context.stroke()

def rect(input, context):
    context.rectangle(input.x+.5, input.y+.5, input.width, input.height)

    if input.color:
        red, green, blue = RGBTuple(input.color)
        context.set_source_rgba(red,green,blue,1)
        context.set_line_width(input.stroke)
        context.stroke_preserve()

    if input.fill:
        r,g,b = RGBTuple(input.fill)
        context.set_source_rgba(r,g,b,1)
        context.fill()


"""
render the conf dist graph to a png and return it as the response
this results in the image being downloaded by the user
"""
def renderToPNG(request):
    import cairo

    if request.method == 'POST': # If the form has been submitted
        form = PlotForm(request.POST) # A form bound to the POST data
        if form.is_valid(): # All validation rules pass
            data = form.cleaned_data
            svg, bins = drawGraph(
                        request,
                        data['x'],
                        data['y'],
                        data['x1'],
                        data['y1'],
                        data['attribute'],
                        data['xProperty'],
                        data['yProperty'],
                        data['reference'],
                        int(data['residue']),
                        data['xBin'],
                        data['yBin'])

    else:
        form = PlotForm() # An unbound form
        svg,bins = drawGraph(request)

    width = 450
    height = 450

    response = HttpResponse(mimetype="image/png")
    response['Content-Disposition'] = 'attachment; filename="plot.png"'
    surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, width, height)
    ctx = cairo.Context (surface)

    for rec in svg.rects:
        rect(rec, ctx)

    for action in svg.lines:
        line(action, ctx)

    for bin in bins:
        svgrec = Rect(bin[0], bin[1], bin[3], bin[2], 1, bin[4], bin[4])
        rect(svgrec, ctx)

    for text in svg.texts:
        ctx.set_source_rgba(0,0,0,1)
        ctx.set_font_size (12);
        ctx.move_to (text.x, text.y);
        ctx.show_text (text.text);


    surface.write_to_png(response)

    return response


"""
render conf dist plot using jquery.svg
"""
def renderToSVG(request):

    response_dict = {
        'referenceValues' : RefDefaults(),
        'stats_fields':STATS_FIELDS
        }

    if request.method == 'POST': # If the form has been submitted
        form = PlotForm(request.POST) # A form bound to the POST data
        if form.is_valid(): # All validation rules pass
            data = form.cleaned_data
            svg, boxes = drawGraph(
                        request,
                        data['x'],
                        data['y'],
                        data['x1'],
                        data['y1'],
                        data['attribute'],
                        data['xProperty'],
                        data['yProperty'],
                        data['reference'],
                        int(data['residue']),
                        data['xBin'],
                        data['yBin'])

            # get values out of the form
            response_dict['xProperty'] = form.cleaned_data['xProperty']
            response_dict['yProperty'] =form.cleaned_data['yProperty']
            response_dict['xBin'] = form.cleaned_data['xBin']
            response_dict['yBin'] = form.cleaned_data['yBin']
            response_dict['attribute'] = form.cleaned_data['attribute']

    else:
        form = PlotForm() # An unbound form
        svg,boxes = drawGraph(request)

        # get default values from the form
        response_dict['xProperty'] = form.fields['xProperty'].initial
        response_dict['yProperty'] = form.fields['yProperty'].initial
        response_dict['xBin'] = form.fields['xBin'].initial
        response_dict['yBin'] = form.fields['yBin'].initial
        response_dict['attribute'] = form.fields['attribute'].initial

    response_dict['form']   = form
    response_dict['svg']    = svg
    response_dict['boxes']  = boxes

    return render_to_response('graph.html', response_dict, context_instance=RequestContext(request))


"""
render the results of the search as a TSV (tab separated file)
and return it to the user as a download
"""
def plotDump(request):
    if request.method == 'POST': # If the form has been submitted
        form = PlotForm(request.POST) # A form bound to the POST data
        if form.is_valid(): # All validation rules pass
            data = form.cleaned_data

            cdp = ConfDistPlot(
                400,            #height
                400,            #width
                0,              #Xpadding
                0,              #Ypadding
                55,              #Xoffset
                45,              #Yoffset
                data['x'],      #Xstart
                data['x1'],           #Xend
                data['y'],         #Ystart
                data['y1'],           #Yend
                data['xBin'],           #Xbin
                data['yBin'],           #Ybin
                data['xProperty'],      #X property
                data['yProperty'],      #Y property
                data['attribute'],#property
                #data['reference'],
                int(data['residue']),
                request.session['search'].querySet()
            )

            response = HttpResponse(mimetype="text/tab-separated-values")
            response['Content-Disposition'] = 'attachment; filename="plot.tsv"'

            # get search out of the session and pass it on
            # TODO ^^

            cdp.Plot()
            cdp.PrintDump(response)

            return response

    return None




