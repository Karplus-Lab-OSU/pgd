from django.conf.urls.defaults import *

urlpatterns = patterns('pgd_search.views',
    (r'^plot/svg/$', 'renderToSVG'),
    (r'^plot/png/$', 'renderToPNG'),
)
