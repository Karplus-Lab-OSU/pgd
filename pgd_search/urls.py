from django.conf.urls import *
from pgd_search.search.views import search, saved, editSearch, help, qtiphelp, saveSearch, deleteSearch, protein_search, chi_help
from pgd_search.plot.views import renderToSVG, renderToPNG, plotDump, plot
from pgd_search.statistics.views import search_statistics, search_statistics_data, search_statistics_aa_data
from pgd_search.dump.views import dataDump
from pgd_search.browse.views import browse
from pgd_search.histogram.views import renderHist

urlpatterns = patterns('',
    (r'^$', search),
    (r'^protein_code/$', protein_search),
    (r'^results/$', plot),
    (r'^plot/svg/$', plot),
    (r'^plot/svg/render/$', renderToSVG),
    (r'^plot/png/$', renderToPNG),
    (r'^plot/dump/$', plotDump),
    (r'^statistics/$', search_statistics),
    (r'^statistics/data/$', search_statistics_data),
    (r'^statistics/aa/(\w+)$', search_statistics_aa_data),
    (r'^dump/$', dataDump),
    (r'^browse/$', browse),
    (r'^histogram/render/$', renderHist),
    (r'^saved/$', saved),
    (r'^help/$', help),
    (r'^chi_help/$', chi_help),
    (r'^qtiphelp/$', qtiphelp),
   
    #old style, for compatibility in the case of forgotton/overlooked links
    (r'^viewSearch/$', editSearch),
    (r'^viewSearch/(?P<search_id>\d+)/$', editSearch),
    (r'^edit/$', editSearch),
    (r'^edit/(?P<search_id>\d+)/$', editSearch),
    (r'^saveSearch/$', saveSearch),
    (r'^saveSearch/(?P<search_id>\d+)/$', saveSearch),
    (r'^deleteSearch/$', deleteSearch),
    (r'^deleteSearch/(?P<search_id>\d+)/$', deleteSearch),

    #new style 
    (r'^view/$', editSearch),
    (r'^(?P<search_id>\d+)/view/$', editSearch),
    (r'^edit/$', editSearch),
    (r'^(?P<search_id>\d+)/$', editSearch),
    (r'^save/$', saveSearch),
    (r'^(?P<search_id>\d+)/save/$', saveSearch),
    (r'^delete/$', deleteSearch),
    (r'^(?P<search_id>\d+)/delete/$', deleteSearch),
)
