# Create your views here.
from registration.backends.default.views import RegistrationView
from forms import UserRegistrationForm as MyCustomRegistrationForm
from models import RegProfile

class MyRegistrationView(RegistrationView):

    form_class= MyCustomRegistrationForm

    def register(self, request, form):
    	new_user = super(MyRegistrationView, self).register(request, form)
        new_profile = RegProfile(user=new_user, 
        	research_summary=form.cleaned_data['research_summary'], 
        	location=form.cleaned_data['location'],
        	timezone=form.cleaned_data['timezone'],
        	country=form.cleaned_data['country'] 
        	)
        new_profile.save()
        return new_user