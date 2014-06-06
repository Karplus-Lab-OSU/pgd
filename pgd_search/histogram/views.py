from django.http import HttpResponse
from django.utils import simplejson

from Histogram import HistogramPlot

def histogram(request, X, Xm, Y, Ym, histoX, histoY, histoZ, histoXr, histoYr, histoZr):
    query = request.session['search'].querySet()
    hp = HistogramPlot(query, X, Xm, Y, Ym, histoX, histoY, histoZ, histoXr, histoYr, histoZr)
    svg = hp.HistoPlot()
    return svg

def renderHist(request):
    data = request.POST
    svg = histogram(
                    request,
                    data["x"],
                    data["x1"],
                    data["y"],
                    data["y1"],
                    data["xRes"],
                    data["yRes"],
                    data["zRes"],
                    data["xResNum"],
                    data["yResNum"],
                    data["zResNum"]
                    )
    json = simplejson.dumps({'svg':svg.to_dict()})
    return HttpResponse(json)
