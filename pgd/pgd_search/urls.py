from django.conf.urls.defaults import *

urlpatterns = patterns('pgd_search.views',
    (r'^$', 'search'),
    (r'^plot/svg/$', 'renderToSVG'),
    (r'^plot/png/$', 'renderToPNG'),
)
