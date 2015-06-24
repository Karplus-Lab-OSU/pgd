from django.conf.urls import include, url, patterns
from django.conf import settings
from django.contrib import admin
from django.contrib.auth import views as auth_views
from django.core.urlresolvers import reverse_lazy
from views import ReferencesView, ContactUsView, NewsView, WelcomeView
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
    #(r'^admin/', include(admin.site.urls)),

    url(r'^search/', include('pgd_search.urls'), name='pgd_search'),

    url(r'^accounts/', include('pgd_core.urls')),
    
    # Static pages:
    (r'^references/$', ReferencesView.as_view()),
    (r'^contactus/$', ContactUsView.as_view()),
    (r'^news/$', NewsView.as_view()),

    #default url
    url(r'^$', WelcomeView.as_view(), name='pgd_home'),
)

#The following is used to serve up local media files like images
#if settings.LOCAL_DEV:
baseurlregex = r'^static/(?P<path>.*)$'
urlpatterns += patterns('',
    (baseurlregex, 'django.views.static.serve', {'document_root':  settings.MEDIA_ROOT}),
    (r'^media/(?P<path>.*)$', 'django.views.static.serve', {'document_root':  settings.MEDIA_ROOT}),
)
