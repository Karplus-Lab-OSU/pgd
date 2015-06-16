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
      
    #override the default urls
    url(r'^accounts/login/$', auth_views.login,
        {'template_name': 'registration/login.html'},
        name='auth_login'),
    
    url(r'^accounts/logout/$',
       auth_views.logout,
       {'template_name' : 'registration/logout.html'},
       name='auth_logout'),

    url(r'^accounts/password/change/$',
       auth_views.password_change,
       {'post_change_redirect': reverse_lazy('auth_password_change_done'),
       'template_name' : 'registration/password_change.html'},
       name='auth_password_change'),

    url(r'^accounts/password/change/done/$',
       auth_views.password_change_done,
       {'template_name' : 'registration/change_done.html'},
       name='auth_password_change_done'),

    url(r'^accounts/password/reset/$',
       auth_views.password_reset,
       {'post_reset_redirect': reverse_lazy('auth_password_reset_done'), 
       'template_name': 'registration/password_reset.html'},
       name='auth_password_reset'),

    url(r'^accounts/password/reset/confirm/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>.+)/$',
       auth_views.password_reset_confirm,
       {'post_reset_redirect': reverse_lazy('auth_password_reset_complete'),
       'template_name' : 'registration/reset_confirm.html'},
       name='auth_password_reset_confirm'),

    url(r'^accounts/password/reset/complete/$',
       auth_views.password_reset_complete,
       {'post_reset_redirect': reverse_lazy('auth_password_reset_complete'), 
        'template_name' : 'reset_complete'},
       name='auth_password_reset_complete'),

    url(r'^accounts/password/reset/done/$',
       auth_views.password_reset_done,
       {'template_name': 'registration/reset_done.html'},
       name='auth_password_reset_done'),

    url(r'^accounts/', include('registration.backends.default.urls')),

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
