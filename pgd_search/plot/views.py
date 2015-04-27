import math
import pickle
from django.db.models import Max, Min
from django.http import HttpResponse
from django.template import RequestContext
from django.conf import settings
from django.shortcuts import render_to_response
import json

from PlotForm import PlotForm, ATTRIBUTE_CHOICES, PROPERTY_CHOICES
from ConfDistFuncs import *
from pgd_constants import AA_CHOICES
from pgd_search.views import settings_processor
from pgd_splicer.sidechain import sidechain_string_dict

AA_CHOICES = [aa[1].upper() for aa in filter(lambda x: x[1].upper() in sidechain_string_dict, AA_CHOICES)]

def drawGraph(request, height=470, width=560, xStart=None, yStart=None, xEnd=None, yEnd=None, attribute='Observations', xProperty='phi', yProperty='psi', reference=None, sigmaVal=3, residue_attribute=None, residue_xproperty=None, residue_yproperty=None, xBin=None, yBin=None, background_color='#ffffff',graph_color='#222222',text_color='#000000', hue='green', hash_color='666666'):
    """
    Renders a conformational distribution graph
    @return: returns an SVG instance.
    """

    query = pickle.loads(request.session['search']).querySet()
    # calculate default values for min, max, and binsize if no values were given
    if residue_xproperty == 0:
        xPrefix = ''
    elif residue_xproperty < 0:
        xPrefix = ''.join(['prev__' for i in range(residue_xproperty, 0)])
    else:
        xPrefix = ''.join(['next__' for i in range(residue_xproperty)])

    if residue_yproperty == 0:
        yPrefix = ''
    elif residue_yproperty < 0:
        yPrefix = ''.join(['prev__' for i in range(residue_yproperty, 0)])
    else:
        yPrefix = ''.join(['next__' for i in range(residue_yproperty)])

    if xStart == None:
        xStart = query.aggregate(min=Min('%s%s' % (xPrefix, xProperty)))['min']
    if xEnd == None:
        xEnd = query.aggregate(max=Max('%s%s' % (xPrefix, xProperty)))['max']
    if yStart == None:
        yStart = query.aggregate(min=Min('%s%s' % (yPrefix ,yProperty)))['min']
    if yEnd == None:
        yEnd = query.aggregate(max=Max('%s%s' % (yPrefix ,yProperty)))['max']
    if xBin == None:
        xBin = math.fabs(xEnd - xStart) / 36
    if yBin == None:
        yBin = math.fabs(yEnd - yStart) / 36

    try:
        cdp = ConfDistPlot(
                width,    #width
                height,   #height
                xStart,         #Xstart
                xEnd,           #Xend
                yStart,         #Ystart
                yEnd,           #Yend
                xBin,           #Xbin
                yBin,           #Ybin
                xProperty,      #X property
                yProperty,      #Y property
                attribute,      #property
                sigmaVal,
                residue_attribute,
                residue_xproperty,
                residue_yproperty,
                pickle.loads(request.session['search']).querySet(),
                hue,
                background_color,
                graph_color,
                text_color,
                hash_color
        )

        svg = cdp.Plot()
    except Exception, e:
        print 'exception', e
        import traceback, sys
        exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        print "*** print_tb:"
        traceback.print_tb(exceptionTraceback, limit=10, file=sys.stdout)

        raise e
    return (svg, xStart, xEnd, xBin, yStart, yEnd, yBin)




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
            svg, x,x1,xBin,y,y1,yBin = drawGraph(
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
                        int(data['sigmaVal']),
                        int(data['residue_attribute']),
                        int(data['residue_xproperty']),
                        int(data['residue_yproperty']),
                        data['xBin'],
                        data['yBin'],
                        data['background_color'],
                        data['graph_color'],
                        data['text_color'],
                        data['plot_hue'],
                        data['hash_color'])

    else:
        form = PlotForm() # An unbound form
        svg,x,x1,xBin,y,y1,yBin = drawGraph(request)
        width = 560
        height = 480

    response = HttpResponse(mimetype="image/png")
    response['Content-Disposition'] = 'attachment; filename="plot.png"'
    svg.render_png(response, width, height+30)

    return response


def plot(request):
    """
    Draws the plot page.  The plot page will rely on AJAX calls to 
    render the graph
    """
    form = PlotForm() # An unbound form
    response_dict = {
        'defaults' : json.dumps(RefDefaults()),
        'xProperty': form.fields['xProperty'].initial,
        'yProperty': form.fields['yProperty'].initial,
        'xBin': form.fields['xBin'].initial,
        'yBin': form.fields['yBin'].initial,
        'attribute': form.fields['attribute'].initial,
        'form': form,
        'attribute_choices':ATTRIBUTE_CHOICES,
        'property_choices':PROPERTY_CHOICES,
        'sidechain_angles':bond_angles_string_dict,
        'sidechain_lengths':bond_lengths_string_dict,
        'aa_choices':AA_CHOICES
    }
    return render_to_response('graph.html', response_dict, context_instance=RequestContext(request, processors=[settings_processor]))


def renderToSVG(request):
    """
    render conf dist plot using jquery.svg
    """
    try:
        form = PlotForm(request.POST) # A form bound to the POST data
        if form.is_valid(): # All validation rules pass
            data = form.cleaned_data
            svg,x,x1,xBin,y,y1,yBin = drawGraph(
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
                                                int(data['sigmaVal']),
                                                int(data['residue_attribute']),
                                                int(data['residue_xproperty']),
                                                int(data['residue_yproperty']),
                                                data['xBin'],
                                                data['yBin'],
                                                data['background_color'],
                                                data['graph_color'],
                                                data['text_color'],
                                                data['plot_hue'],
                                                data['hash_color'])
            _json = json.dumps({'svg':svg.to_dict(), \
                                        'x':x, 'x1':x1, 'xBin':xBin, \
                                        'y':y, 'y1':y1, 'yBin':yBin})
            return HttpResponse(_json)

        else:
            """
            Errors in the form - repackage the error list as a list of errors
            This list can then be json serialized and processed by the javascript
            on the plot page
            """
            errors = []
            for k, v in form.errors.items():
                for error in v:
                    errors.append([k, error._proxy____args[0]])

            return HttpResponse(json.dumps({'errors':errors}))
    except Exception, e:
        print 'exception', e
        import traceback, sys
        exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        print "*** print_tb:"
        traceback.print_tb(exceptionTraceback, limit=10, file=sys.stdout)
        return HttpResponse("-1")


def plotDump(request):
    """
    render the results of the search as a TSV (tab separated file)
    and return it to the user as a download
    """
    if request.method == 'POST': # If the form has been submitted
        form = PlotForm(request.POST) # A form bound to the POST data
        if form.is_valid(): # All validation rules pass
            data = form.cleaned_data

            cdp = ConfDistPlot(
                360,               #height
                360,               #width
                data['x'],         #Xstart
                data['x1'],        #Xend
                data['y'],         #Ystart
                data['y1'],        #Yend
                data['xBin'],      #Xbin
                data['yBin'],      #Ybin
                data['xProperty'], #X property
                data['yProperty'], #Y property
                #data['attribute'],#property
                'all',#property
                #data['reference'],
                int(data['sigmaVal']),
                int(data['residue_attribute']),
                int(data['residue_xproperty']),
                int(data['residue_yproperty']),
                pickle.loads(request.session['search']).querySet()
            )

            response = HttpResponse(mimetype="text/tab-separated-values")
            response['Content-Disposition'] = 'attachment; filename="plot.tsv"'

            cdp.Plot()
            cdp.PrintDump(response)

            return response

    return HttpResponse('Error')



