from django.shortcuts import render_to_response
from django.template import RequestContext
from pgd_search.views import settings_processor
from pgd_splicer.models import pdb_select_settings
import settings

extra_context = {
    'SITE_ROOT': settings.SITE_ROOT,
    'PGD_VERSION': settings.PGD_VERSION,
    'ROOT': settings.SITE_ROOT,
    'DATA_VERSION': pdb_select_settings.DATA_VERSION,
    'GOOGLE_ID': settings.GOOGLE_ID
}


def welcome(request):
	#delete session variable
	try:
		del request.session['search']
	except:
		pass
	return render_to_response('welcome.html',
						{'extra_context' : extra_context},
						context_instance=RequestContext(request, processors=[settings_processor]))

def references(request):
	#delete session variable
	try:
		del request.session['search']
	except:
		pass
	return render_to_response('references.html',
						{'extra_context' : extra_context},
						context_instance=RequestContext(request, processors=[settings_processor]))

def news(request):
	#delete session variable
	try:
		del request.session['search']
	except:
		pass
	return render_to_response('news.html',
						{'extra_context' : extra_context},
						context_instance=RequestContext(request, processors=[settings_processor]))

def contactus(request):
	#delete session variable
	try:
		del request.session['search']
	except:
		pass
	return render_to_response('contactus.html',
						{'extra_context' : extra_context},
						context_instance=RequestContext(request, processors=[settings_processor]))
