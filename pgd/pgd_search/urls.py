from django.conf.urls.defaults import *

<<<<<<< HEAD:pgd/pgd_search/urls.py
from pgd_search.search.views import search, saved, editSearch, help, qtiphelp
from pgd_search.plot.views import renderToSVG, renderToPNG, plotDump
=======
from pgd_search.search.views import search, saved, editSearch, help
from pgd_search.plot.views import renderToSVG, plot, renderToPNG, plotDump
>>>>>>> 341a8ed50b4c0f217315479f764bd4d4898d29c2:pgd/pgd_search/urls.py
from pgd_search.statistics.views import searchStatistics
from pgd_search.dump.views import dataDump
from pgd_search.browse.views import browse

urlpatterns = patterns('',
    (r'^$', search),
    (r'^results/$', plot),
    (r'^plot/svg/$', plot),
    (r'^plot/svg/render/$', renderToSVG),
    (r'^plot/png/$', renderToPNG),
    (r'^plot/dump/$', plotDump),
    (r'^statistics/$', searchStatistics),
    (r'^dump/$', dataDump),
    (r'^browse/$', browse),
    (r'^saved/$', saved),
    (r'^help/$', help),
	(r'^qtiphelp/$', qtiphelp),
    (r'^edit/$', editSearch),
    (r'^edit/(?P<search_id>\d+)/$', editSearch),
)
