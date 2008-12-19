from django.conf.urls.defaults import *

from pgd_search.search.views import search
from pgd_search.plot.views import renderToSVG, renderToPNG, plotDump
from pgd_search.statistics.views import searchStatistics
from pgd_search.dump.views import dataDump
from pgd_search.browse.views import browse

urlpatterns = patterns('',
    (r'^$', search),
    (r'^results/$', renderToSVG),
    (r'^plot/svg/$', renderToSVG),
    (r'^plot/png/$', renderToPNG),
    (r'^plot/dump/$', plotDump),
    (r'^statistics/$', searchStatistics),
    (r'^dump/$', dataDump),
    (r'^browse/$', browse),
)
