"""
This file demonstrates writing tests using the unittest module. These will pass
when you run "manage.py test".

Replace this with more appropriate tests for your application.
"""

from django.test import TestCase , Client
from django.core.urlresolvers import reverse


class SimpleTest(TestCase):


	fixtures=['users.json']

	def test_reset_password(self):
		respon = self.client.post('http://testserver/accounts/password/reset/', {'email' : 'email@gmail.com'})
		self.assertEqual(respon.status_code, 302)	
		self.assertEqual(respon['Location'], 'http://testserver/accounts/password/reset/done/')

	def test_change_password(self):
		resp1 = Client()
		resp1_get = resp1.get('http://testserver/accounts/password/change/')
		self.assertEqual(resp1_get.status_code, 302)
		resp2_post1 = resp1.post(resp1_get['Location'], {'username':'test1', 'password':'hey'}, follow=True)
		self.assertEqual(resp2_post1.status_code, 200)
		resp2_post2 = resp1.post(resp2_post1.redirect_chain[0][0], {'old_password': 'hey',
		'new_password1': 'heythere', 'new_password2' : 'heythere'}, follow=True)
		self.assertEqual(resp2_post2.status_code, 200)
		self.assertEqual(resp2_post2.redirect_chain[0][0], 'http://testserver/accounts/password/change/done/')	