# Create your views here.
from django.shortcuts import render
from django.shortcuts import redirect
from django.core.urlresolvers import reverse
from django.contrib.auth.models import User
from django.contrib.auth.decorators import login_required
from registration.backends.default.views import RegistrationView
from forms import UserRegistrationForm as MyCustomRegistrationForm, EditForm

class MyRegistrationView(RegistrationView):

    form_class= MyCustomRegistrationForm

    def register(self, request, form):
    	new_user = super(MyRegistrationView, self).register(request, form)
        new_user.first_name = form.cleaned_data['first_name']
        new_user.last_name  = form.cleaned_data['last_name']
        new_user.save()
        return new_user


@login_required(login_url='/accounts/login')
def profile_view(request) :
	
	if request.user.is_active :
		#Fetch information for this user from RegProfile
		prof_details = {}

		prof_details['first_name'] = request.user.first_name
		prof_details['last_name'] = request.user.last_name
		prof_details['email'] = request.user.email
		prof_details['user_name'] = request.user.username
		prof_details['full_name'] = request.user.get_full_name()

		return render(request, 'profile.html', prof_details)
	else :
		return redirect('/accounts/login')


@login_required(login_url='/accounts/login')
def edit_profile_view(request):

	if request.user.is_active :
		if request.method == 'POST':
			form = EditForm(request.POST)
			if form.is_valid():
				user = request.user
				user.first_name = form.cleaned_data['first_name']
				user.last_name  = form.cleaned_data['last_name']
				user.save()
				return render(request, 'profile_edited.html')
		else:
			form = EditForm()
	else :
		return redirect('/accounts/login')
	return render(request, 'edit_profile.html', {'form' : form,}) 


def get_profile_view(request, username):

	try:
		user = User.objects.get(username=username)
		prof_details = {}
		prof_details['first_name'] = user.first_name
		prof_details['last_name'] = user.last_name
		prof_details['email'] = user.email
		prof_details['user_name'] = username
		if request.user.is_active :
			prof_details['full_name'] = request.user.get_full_name()
		else :
			prof_details['full_name'] = ''
		return render(request, 'profile.html', prof_details)
	
	except Exception, e:
		return render(request, 'user_non-existant.html')	

		#raise e
	