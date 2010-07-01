import math

from django.db.models import Max, Min
from django.http import HttpResponse
from django.template import RequestContext
from django.conf import settings
from django.shortcuts import render_to_response
from django.utils import simplejson

from Histogram import *
from pgd_search.views import settings_processor

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