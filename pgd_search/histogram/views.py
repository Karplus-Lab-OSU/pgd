from django.http import HttpResponse
import json
import pickle
from Histogram import HistogramPlot

def histogram(request, X, Xm, Y, Ym, histoX, histoY, histoZ, histoXr, histoYr, histoZr):
    query = pickle.loads(request.session['search']).querySet()
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
    _json = json.dumps({'svg':svg.to_dict()})
    return HttpResponse(_json)
