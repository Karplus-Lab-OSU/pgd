from django.http import HttpResponse
from django.conf import settings
from django.shortcuts import render_to_response

from PlotForm import PlotForm
from ConfDistFuncs import *
from svg import *



"""
Renders a conformational distribution graph
@return: retusns an SVG instance.
"""
def drawGraph(xStart=-180, yStart=-180, xEnd=180, yEnd=180, attribute='Observations', xProperty='phi', yProperty='psi', reference=None, residue=None, xBin=10, yBin=10):
    svg = SVG()

    x = 55;
    y = 45;
    height = 400;
    width = 400;
    hashsize = 10

    #background
    svg.rect(x, y, width, height, 1, '#00fffff', '#222222');
    #svg.rect(0, 0, width+90, height+90, 1, '#00fffff');
    #svg.rect(x, y, width, height, 1, '#666666');

    #border
    svg.rect(x, y, width, height, 1, '#000000');

    #axis
    svg.line( x, y+height/2, x+width, y+height/2, 1, '#666666');
    svg.line( x+width/2, y, x+width/2, y+height, 1, '#666666');

    #hashes
    for i in range(9):
        hashx = x+(width/8.0)*i
        hashy = y+(height/8.0)*i
        svg.line( hashx, y+height, hashx, y+height+hashsize, 1, '#000000');
        svg.line( x, hashy, x-hashsize, hashy, 1, '#000000');

    #x axis text
    xtext  = xStart
    xtext1 = xEnd
    step = (xtext1 - xtext) / 4
    for i in range(5):
        text = xtext + step*i
        hashx = x+(width/4)*i-(2.5*len(str(text)))
        svg.text(hashx, y+height+hashsize*2+3, str(text),12)

    #y axis text
    ytext  = yStart
    ytext1 = yEnd
    step = (ytext1 - ytext) / 4
    for i in range(5):
        text = ytext + step*i
        hashy = y+(height/4)*i+4
        svg.text(x-5-(8*len(str(text))), hashy, str(text),12)

    #title text
    len1 = 220 - len(xProperty)*7/2 - len(xProperty)*7/2
    len2 = 182 - len(attribute)*7/2
    svg.text(len1,15, 'Plot of %s vs. %s' % (xProperty,yProperty), 12)
    svg.text(len2,35, 'Shading Based Off of %s' % attribute, 12)

    cdp = ConfDistPlot(
            400,            #height
            400,            #width
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
            residue         #residue Index
            #reference
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
        svg,bins = drawGraph()

    width = 500
    height = 500

    response = HttpResponse(mimetype="image/png")
    response['Content-Disposition'] = 'attachment; filename="plot.png"'
    surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, width, height)
    ctx = cairo.Context (surface)

    for rec in svg.rects:
        rect(rec, ctx)

    for action in svg.lines:
        line(action, ctx)

    for bin in bins:
        svgrec = Rect(bin[0], bin[1], bin[2], bin[3], 1, bin[4], bin[4])
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
    if request.method == 'POST': # If the form has been submitted
        form = PlotForm(request.POST) # A form bound to the POST data
        if form.is_valid(): # All validation rules pass
            data = form.cleaned_data
            svg, boxes = drawGraph(
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
        svg,boxes = drawGraph()

    return render_to_response('graph.html', {
        'SITE_ROOT': settings.SITE_ROOT,
        'MEDIA_URL': settings.MEDIA_URL,
        'form': form,
        'svg': svg,
        'boxes': boxes,
        'referenceValues' : RefDefaults()
    })


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
                int(data['residue'])
            )

            response = HttpResponse(mimetype="text/tab-separated-values")
            response['Content-Disposition'] = 'attachment; filename="plot.tsv"'

            # get search out of the session and pass it on
            # TODO ^^

            cdp.Plot()
            cdp.PrintDump(response)

            return response

    return None




