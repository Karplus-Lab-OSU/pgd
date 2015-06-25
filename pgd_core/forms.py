import pytz
from django import forms
from registration.forms import RegistrationForm
from django.utils.translation import ugettext_lazy as _
from pgd_core.models import RegProfile
from registration.models import RegistrationProfile

attrs_dict = { 'class': 'required' }

class UserRegistrationForm(RegistrationForm):


	TIMEZONE = pytz.all_timezones


	research_summary = forms.CharField(
						widget=forms.Textarea(attrs={'class': 'research'}))

	location         = forms.CharField(label="university/company", 
						widget=forms.TextInput(attrs={'placeholder': 'Search'}))

	timezone 	 	 = forms.CharField(
						widget=forms.TextInput(attrs={'class': 'timezone'}) )

	country			 = forms.CharField(widget=forms.TextInput(attrs={'class': 'country'}))

'''
	def save(self, profile_callback=None):
		new_user = RegistrationProfile.objects.create_inactive_user(username=self.cleaned_data['username'],
			password=self.cleaned_data['password1'],
			email=self.cleaned_data['email'])
		new_profile = UserProfile(user=new_user, research_summary=self.cleaned_data['research_summary'], 
			location=self.cleaned_data['location'], timezone=self.cleaned_data['timezone'], 
			country=self.cleaned_data['country'])
		new_profile.save()
		return new_user
'''