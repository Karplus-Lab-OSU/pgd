from django.http import HttpResponse
from django.template import RequestContext, Context, loader
from django.conf import settings

from ConfDistFuncs import *
from svg import *
#from pgd_core.models import *

""" 
Renders a conformational distribution graph
@return: retusns an SVG instance.
"""
def drawGraph():
    svg = SVG()

    x = 100;
    y = 100;
    height = 400;
    width = 400;
    hashsize = 10

    #background
    svg.rect(x, y, width, height, 1, '#00fffff', '#222222');
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
    xtext  = -180
    xtext1 = 180
    step = (xtext1 - xtext) / 4
    for i in range(5):
        text = xtext + step*i
        hashx = x+(width/4)*i-(4*len(str(text)))
        svg.text(hashx,y+height+hashsize*3,str(text),14)

    #y axis text
    ytext  = -180
    ytext1 = 180
    step = (ytext1 - ytext) / 4
    for i in range(5):
        text = xtext + step*i
        hashy = y+(height/4)*i+7
        svg.text(x-20-(8*len(str(text))),hashy,str(text),14)

    cdp = ConfDistPlot(
            400,            #height
            400,            #width
            0,             #Xpadding
            0,             #Ypadding
            100,              #Xoffset
            100,              #Yoffset
            -180,           #Xstart
            180,            #Xend
            -180,           #Ystart
            180,            #Yend
            10,             #Xbin
            10,             #Ybin
            "phi",          #X property
            "psi",          #Y property
            '1sny',         #protein code
            'Observations'  #property
    )

    svg = cdp.Plot(svg)
    return svg

def RGBTuple(rgbString):
    sub = rgbString[-6:]
    red = int(sub[:2],16)/255.0
    green = int(sub[2:4],16)/255.0
    blue = int(sub[4:], 16)/255.0
    return (red,green,blue)

def line(input, context):
    context.move_to(input.x, input.y)
    context.line_to(input.x1, input.y1)
    context.set_line_width(input.stroke)
    r,g,b = RGBTuple(input.color)
    context.set_source_rgba(r,g,b,1)
    context.stroke()

def rect(input, context):
    context.rectangle(input.x, input.y, input.width, input.height)

    if input.fill:
        r,g,b = RGBTuple(input.fill)
        context.set_source_rgba(r,g,b,1)
        context.fill()

    if input.color:
        red, green, blue = RGBTuple(input.color)
        context.set_source_rgba(red,green,blue,1)
        context.set_line_width(input.stroke)
        context.stroke()

"""
render the conf dist graph to a png and return it as the response
this results in the image being downloaded by the user
"""
def renderToPNG(request):
    import cairo

    width = 600
    height = 600

    response = HttpResponse(mimetype="image/png")

    surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, width, height)
    ctx = cairo.Context (surface)

    svg = drawGraph()
    for rec in svg.rects:
        rect(rec, ctx)

    for action in svg.lines:
        line(action, ctx)

    surface.write_to_png(response)

    return response

"""
render conf dist plot using jquery.svg
"""
def renderToSVG(request):
    svg = drawGraph()

    t = loader.get_template('graph.html')
    c = RequestContext(request, {
        'MEDIA_URL': settings.MEDIA_URL,
        'svg': svg
    })

    return HttpResponse(t.render(c))