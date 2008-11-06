from django.conf.urls.defaults import *

urlpatterns = patterns('pgd.tasks.views',
    (r'^$', 'showtasks'),
    (r'^progress/$', 'taskprogress'),
    (r'^start/$', 'starttask'),
)
