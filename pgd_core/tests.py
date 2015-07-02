"""
This file demonstrates writing tests using the unittest module. These will pass
when you run "manage.py test".

Replace this with more appropriate tests for your application.
"""

from django.test import TestCase , Client
from django.core.urlresolvers import reverse
from django.core import mail


class RegistrationTestCase(TestCase):


	fixtures    	  = ['pgd_core']
	default_url 	  = "http://testserver"
	test_email		  = {'email' : 'email@example.org'}
	test_credentials = {'username':'test_user', 'password':'hello'}
	
	register_details  = {
	'username' : 'vamos', 
	'email' : 'email@example.org', 
	'password1' : 'hola',
	'password2' : 'hola',
	'first_name' : 'Test',
	'last_name' : 'User'}

	test_new_password = {'old_password': 'hello','new_password1': 'heythere',
						'new_password2' : 'heythere'}

	change_name = {'first_name' : 'Dummy',
					'last_name' : 'User'}

	def test_register(self) :
		test_client = Client()
		response = test_client.post(reverse('registration_register'), self.register_details, follow=True)
		self.assertEqual(response.status_code, 200)
		self.assertEqual(response.redirect_chain[-1][0], self.default_url+'/accounts/register/complete/')
		self.assertEqual(len(mail.outbox), 1)
		self.assertEqual(mail.outbox[0].to[0], self.test_email['email'])	

	def test_reset_password(self):
		response = self.client.post(reverse('auth_password_reset'), self.test_email)
		self.assertEqual(response.status_code, 302)
		self.assertEqual(response['Location'], self.default_url+reverse('auth_password_reset_done'))


	def test_edit_user_profile(self):

		test_client = Client()
		get_response = test_client.get(reverse('user_profile_edit'), follow=True)
		self.assertEqual(get_response.status_code, 200)
		post_response_1 = test_client.post(get_response.redirect_chain[-1][0], self.test_credentials, follow=True)
		self.assertEqual(post_response_1.status_code, 200)
		edit_post = test_client.post(post_response_1.redirect_chain[-1][0], self.change_name, follow=True)
		from django.contrib.auth.models import User
		user = User.objects.get(username='test_user')
		self.assertEqual(user.first_name, self.change_name['first_name'])
		self.assertEqual(user.last_name, self.change_name['last_name'])



	def test_user_profile(self) :

		#Try to access the user Profile first, before login
		test_client = Client()
		get_profile = test_client.get(reverse('user_profile'), follow=True)
		self.assertEqual(get_profile.status_code, 200)
		post_credentials = test_client.post(get_profile.redirect_chain[-1][0], self.test_credentials, follow=True)
		self.assertEqual(post_credentials.status_code, 200)
		self.assertEqual(post_credentials.redirect_chain[-1][0], self.default_url+reverse('user_profile'))


	def test_change_password(self):

		test_client = Client()
		get_response = test_client.get(reverse('auth_password_change'))
		self.assertEqual(get_response.status_code, 302)
		post_response_1 = test_client.post(get_response['Location'], self.test_credentials, follow=True)
		self.assertEqual(post_response_1.status_code, 200)
		post_response_2 = test_client.post(post_response_1.redirect_chain[-1][0], self.test_new_password, follow=True)
		self.assertEqual(post_response_2.status_code, 200)
		self.assertEqual(post_response_2.redirect_chain[-1][0], self.default_url+reverse('auth_password_change_done'))

	def test_save_search(self):

		test_client = Client()
		get_response = test_client.get(reverse('generic_profile', args=('test_user',)))
		self.assertNotIn('<td class="sSearch">False</td>' , get_response)
		get_profile = test_client.get(reverse('user_profile'), follow=True)
		self.assertEqual(get_profile.status_code, 200)
		post_credentials = test_client.post(get_profile.redirect_chain[-1][0], self.test_credentials, follow=True)
		self.assertEqual(post_credentials.status_code, 200)
		profile_page = test_client.get(reverse('generic_profile', args=('test_user',)))
		self.assertIn('<td class="sSearch">False</td>' , profile_page.content)

	def test_search_user(self):

		test_client = Client()
		#Multiple match
		search_multiple = test_client.get(reverse('user-search'), {'q':'test'})
		self.assertIn('<a href="/accounts/profile/test/"> test </a>' , search_multiple.content)
		self.assertIn('a href="/accounts/profile/test_user/"> test_user </a>' , search_multiple.content)

		#Single match 
		search_single = test_client.get(reverse('user-search'), {'q' : 'test_user'})
		self.assertIn('a href="/accounts/profile/test_user/"> test_user </a>' , search_single.content)

		#No matches
		search_none = test_client.get(reverse('user-search'), {'q' : 'whatever'}, redirect=True)
		self.assertEqual(search_none['Location'], self.default_url+reverse('notfound'))
