from django.conf.urls.defaults import *
from django.views.generic.simple import direct_to_template
import settings
from django.contrib import admin
admin.autodiscover()

from pgd import VERSION
from pgd_splicer.models import pdb_select_settings

extra_context = {'SITE_ROOT':settings.SITE_ROOT,
                'MEDIA_ROOT': settings.MEDIA_URL,
                'PGD_VERSION':VERSION,
                'MEDIA':settings.MEDIA_URL,
                'ROOT':settings.SITE_ROOT,
                'DATA_VERSION':pdb_select_settings.DATA_VERSION,
                'GOOGLE_ID':settings.GOOGLE_ID
                }

urlpatterns = patterns('',
    # Example:
    # (r'^pgd/', include('pgd.foo.urls')),

    # Uncomment the admin/doc line below and add 'django.contrib.admindocs' 
    # to INSTALLED_APPS to enable admin documentation:
    # (r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    (r'^admin/(.*)', admin.site.root),
    (r'^accounts/', include('registration_local.urls')),
    (r'^search/', include('pgd_search.urls')),

    # Static pages:
    (r'^references/$','django.views.generic.simple.direct_to_template',
                        {'template':'references.html',
                         'extra_context': extra_context}),
    (r'^contactus/$','django.views.generic.simple.direct_to_template',
                        {'template':'contactus.html',
                         'extra_context': extra_context}),
    (r'^news/$','django.views.generic.simple.direct_to_template',
                        {'template':'news.html',
                         'extra_context': extra_context}),
    #default url
    (r'^$','django.views.generic.simple.direct_to_template',
                        {'template':'welcome.html',
                         'extra_context': extra_context}),
)

#The following is used to serve up local media files like images
#if settings.LOCAL_DEV:
baseurlregex = r'^static/(?P<path>.*)$'
urlpatterns += patterns('',
    (baseurlregex, 'django.views.static.serve', {'document_root':  settings.MEDIA_ROOT}),
)
