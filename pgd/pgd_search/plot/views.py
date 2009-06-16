import cairo
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
def drawGraph(request, height=470, width=470, xStart=-180.0, yStart=-180.0, xEnd=180.0, yEnd=180.0, attribute='Observations', xProperty='phi', yProperty='psi', reference=None, residue=None, xBin=10, yBin=10, background_color='#ffffff',graph_color='#222222',text_color='#000000', hue='green', hash_color='#666666'):

    #store functions locally for speed optimization
    local_len = len
    local_int = int

    svg = SVG()

    #size ratio (470 = 1)
    ratio = width/470.0

    x = round(width*.117);
    y = round(width*.117);
    graph_height = height-2*y;
    graph_width = width-2*x;
    hashsize = 10*ratio

    #image background
    svg.rect(0, 0, width, height, 1, background_color, background_color);

    #graph background
    svg.rect(x, y, graph_width, graph_height, 1, hash_color, graph_color);

    #border
    svg.rect(x, y, graph_width, graph_height, 1, text_color);

    #axis
    if xStart < 0 and xEnd > 0:
        xZero = (graph_width/(xEnd-xStart)) * abs (xStart)
        svg.line( x+xZero, y, x+xZero, y+graph_height, 1, hash_color);

    if yStart < 0 and xEnd > 0:
        yZero = graph_height+y - (graph_height/(yEnd-yStart)) * abs (yStart)
        svg.line( x, yZero, x+graph_width, yZero, 1, hash_color);

    #hashes
    for i in range(9):
        hashx = x+(graph_width/8.0)*i
        hashy = y+(graph_height/8.0)*i
        svg.line( hashx, y+graph_height, hashx, y+graph_height+hashsize, 1, hash_color);
        svg.line( x, hashy, x-hashsize, hashy, 1, hash_color);

    #create a cairo surface to calculate text sizes
    surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, width, height)
    ctx = cairo.Context (surface)
    ctx.set_font_size (12*ratio);

    #labels
    xstep = ((xEnd - xStart)%360 if xProperty in ANGLES else (xEnd - xStart))/ 4
    if not xstep: xstep = 90
    ystep = ((yEnd - yStart)%360 if yProperty in ANGLES else (yEnd - yStart))/ 4
    if not ystep: ystep = 90

    #get Y coordinate for xaxis hashes, this is the same for all x-labels
    xlabel_y = y+graph_height+hashsize*2+(3*ratio)
    for i in range(5):
        #text value
        xtext = ((xStart + xstep*i + 180)%360 - 180) if xProperty in ANGLES else (xStart + xstep*i + 180)
        #drop decimal if value is an integer
        xtext = '%i' % local_int(xtext) if not xtext%1 else '%.1f' %  xtext
        #get X coordinate of hash, offsetting for length of text
        xbearing, ybearing, twidth, theight, xadvance, yadvance = ctx.text_extents(xtext)
        xlabel_x = x+(graph_width/4)*i-xbearing-twidth/2
        #create label
        svg.text(xlabel_x, xlabel_y, xtext,12*ratio, text_color)

        #text value
        ytext = ((yStart + ystep*i + 180)%360 - 180) if yProperty in ANGLES else (yStart + ystep*i + 180)
        #drop decimal if value is an integer
        ytext = '%i' % local_int(ytext) if not ytext%1 else '%.1f' % ytext
        #get Y coordinate offsetting for height of text
        xbearing, ybearing, twidth, theight, xadvance, yadvance = ctx.text_extents(ytext)
        ylabel_y = y+(graph_height/4)*i+(4*ratio)
        #Get X coordinate offsetting for length of hash and length of text
        ylabel_x = (x-(ratio*15))-xbearing-twidth
        #create label
        svg.text(ylabel_x, ylabel_y, ytext,12*ratio, text_color)

    #title text
    title = 'Plot of %s vs. %s' % (xProperty,yProperty)
    xbearing, ybearing, twidth, theight, xadvance, yadvance = ctx.text_extents(title)
    title_x = (width/2) - xbearing - twidth/2
    svg.text(title_x,15*ratio, title, 12*ratio, text_color)

    title = 'Shading Based Off of %s' % attribute
    xbearing, ybearing, twidth, theight, xadvance, yadvance = ctx.text_extents(title)
    title_x = (width/2) - xbearing - twidth/2
    svg.text(title_x,35*ratio, title, 12*ratio, text_color)

    cdp = ConfDistPlot(
            graph_width,    #width
            graph_height,   #height
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
            residue,        #residue Index
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

    if input.color and input.color <> 'None':
        red, green, blue = RGBTuple(input.color)
        context.set_source_rgba(red,green,blue,1)
        context.set_line_width(input.stroke)
        context.stroke_preserve()

    if input.fill and input.fill <> 'None':
        r,g,b = RGBTuple(input.fill)
        context.set_source_rgba(r,g,b,1)
        context.fill()


"""
render the conf dist graph to a png and return it as the response
this results in the image being downloaded by the user
"""
def renderToPNG(request):
    if request.method == 'POST': # If the form has been submitted
        form = PlotForm(request.POST) # A form bound to the POST data
        if form.is_valid(): # All validation rules pass
            data = form.cleaned_data
            width = data['width']
            height = data['height']
            svg, bins = drawGraph(
                        request,
                        height,
                        width,
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
                        data['yBin'],
                        data['background_color'],
                        data['graph_color'],
                        data['text_color'],
                        data['plot_hue'],
                        data['hash_color'])

    else:
        form = PlotForm() # An unbound form
        svg,bins = drawGraph(request)
        width = 460
        height = 460

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
        red, green, blue = RGBTuple(text.color)
        ctx.set_source_rgba(red,green,blue,1)
        ctx.set_font_size (text.size);
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
        }

    if request.method == 'POST': # If the form has been submitted
        form = PlotForm(request.POST) # A form bound to the POST data
        if form.is_valid(): # All validation rules pass
            data = form.cleaned_data
            svg, boxes = drawGraph(
                        request,
                        int(data['height']),
                        int(data['width']),
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
                        data['yBin'],
                        data['background_color'],
                        data['graph_color'],
                        data['text_color'],
                        data['plot_hue'],
                        data['hash_color'])

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

    response_dict['form']         = form
    response_dict['svg']          = svg
    response_dict['boxes']        = boxes

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
                #data['attribute'],#property
                'all',#property
                #data['reference'],
                int(data['residue']),
                request.session['search'].querySet()
            )

            response = HttpResponse(mimetype="text/tab-separated-values")
            response['Content-Disposition'] = 'attachment; filename="plot.tsv"'

            cdp.Plot()
            cdp.PrintDump(response)

            return response

    return None




