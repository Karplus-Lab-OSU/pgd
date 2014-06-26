from django.conf.urls import *
#from django.views.generic.simple import direct_to_template
import settings
from django.contrib import admin
admin.autodiscover()

#from pgd import VERSION
#from pgd_splicer.models import pdb_select_settings

urlpatterns = patterns('',
    # Example:
    # (r'^pgd/', include('pgd.foo.urls')),

    # Uncomment the admin/doc line below and add 'django.contrib.admindocs' 
    # to INSTALLED_APPS to enable admin documentation:
    # (r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    #(r'^admin/(.*)', admin.site.root),
    (r'^accounts/', include('registration.backends.default.urls')),
    (r'^search/', include('pgd_search.urls')),

    # Static pages:
    (r'^references/$', 'pgd.views.references'),
    (r'^contactus/$', 'pgd.views.contactus'),
    (r'^news/$', 'pgd.views.news'),
    
    #default url
    (r'^$','pgd.views.welcome'),
)

#The following is used to serve up local media files like images
#if settings.LOCAL_DEV:
baseurlregex = r'^static/(?P<path>.*)$'
urlpatterns += patterns('',
    (baseurlregex, 'django.views.static.serve', {'document_root':  settings.MEDIA_ROOT}),
    (r'^media/(?P<path>.*)$', 'django.views.static.serve', {'document_root':  settings.MEDIA_ROOT}),
)
