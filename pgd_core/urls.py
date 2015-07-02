from django.conf.urls import *
from django.conf import settings
from django.contrib.auth import views as auth_views
from django.core.urlresolvers import reverse_lazy
from forms import UserRegistrationForm  
from views import MyRegistrationView , profile_view, edit_profile_view, get_profile_view, search, notfound
from registration.backends.default.views import RegistrationView 

urlpatterns = patterns('',

    url(r'^logout/$',
       auth_views.logout,
       {'template_name' : 'registration/logout.html'},
       name='auth_logout'),

    url(r'^password/change/$',
       auth_views.password_change,
       {'post_change_redirect': reverse_lazy('auth_password_change_done'),
       'template_name' : 'registration/password_change.html'},
       name='auth_password_change'),

    url(r'^password/change/done/$',
       auth_views.password_change_done,
       {'template_name' : 'registration/change_done.html'},
       name='auth_password_change_done'),

    url(r'^password/reset/$',
       auth_views.password_reset,
       {'post_reset_redirect': reverse_lazy('auth_password_reset_done'), 
       'template_name': 'registration/password_reset.html'},
       name='auth_password_reset'),

    url(r'^password/reset/confirm/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>.+)/$',
       auth_views.password_reset_confirm,
       {'post_reset_redirect': reverse_lazy('auth_password_reset_complete'),
       'template_name' : 'registration/reset_confirm.html'},
       name='auth_password_reset_confirm'),

    url(r'^password/reset/complete/$',
       auth_views.password_reset_complete,
       {'post_reset_redirect': reverse_lazy('auth_password_reset_complete'), 
        'template_name' : 'reset_complete'},
       name='auth_password_reset_complete'),

    url(r'^password/reset/done/$',
       auth_views.password_reset_done,
       {'template_name': 'registration/reset_done.html'},
       name='auth_password_reset_done'),
    
    url(r'^register/$', MyRegistrationView.as_view(),
        name='registration_register'),

    (r'^', include('registration.backends.default.urls')),
    url(r'^profile/$',  profile_view, name='user_profile'),
    url(r'^profile/([a-zA-Z_@\+\.-]+)/$',get_profile_view, name='generic_profile'),
    url(r'^profile-edit/$',  edit_profile_view, name='user_profile_edit'),
    url(r'^search/$', search, name='user-search'),
    url(r'^notfound/$', notfound ,name='notfound' ),
    )