from django.conf.urls.defaults import *
from django.views.generic.simple import direct_to_template
from registration.views import activate
from registration.views import register
from django.contrib.auth import authenticate, logout
import settings
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    # Example:
    # (r'^pgd/', include('pgd.foo.urls')),

    # Uncomment the admin/doc line below and add 'django.contrib.admindocs' 
    # to INSTALLED_APPS to enable admin documentation:
    # (r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    (r'^admin/(.*)', admin.site.root),
    (r'^settings/', include('dbsettings.urls')),
	(r'^accounts/', include('registration.urls')),
    (r'^search/', include('pgd_search.urls')),

    #default url
    (r'^$','django.views.generic.simple.direct_to_template',
                        {'template':'welcome.html',
                         'extra_context': {'SITE_ROOT':settings.SITE_ROOT,
                                           'MEDIA_ROOT': settings.MEDIA_URL}}),
)

'''urlpatterns += patterns('',
  (r'^accounts/profile/$', direct_to_template, {'template': 'registration/profile.html'}),
  (r'^accounts/password_reset/$', password_reset, {'template_name': 'registration/password_reset.html'}),
  (r'^accounts/password_reset_done/$', password_reset_done, {'template_name': 'registration/password_reset_done.html'}),
  (r'^accounts/password_change/$', password_change, {'template_name': 'registration/password_change.html'}),
  (r'^accounts/password_change_done/$', password_change_done, {'template_name': 'registration/password_change_done.html'}),
)'''

#The following is used to serve up local media files like images
#if settings.LOCAL_DEV:
baseurlregex = r'^static/(?P<path>.*)$'
urlpatterns += patterns('',
    (baseurlregex, 'django.views.static.serve', {'document_root':  settings.MEDIA_ROOT}),
)
