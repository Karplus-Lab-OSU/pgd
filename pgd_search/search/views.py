import math
import re
import pickle

from django.http import HttpResponse, HttpResponseRedirect
from django.template import RequestContext
from django.shortcuts import render_to_response, redirect
from django.conf import settings
from django.forms.util import ErrorList
from django.core.paginator import Paginator, InvalidPage, EmptyPage
import json
from datetime import datetime
from pgd_core.models import Protein
from pgd_search.models import Search, Search_code
from pgd_search.views import RESIDUE_INDEXES, settings_processor
from SearchForm import SearchSyntaxField, SearchForm
from pgd_constants import AA_CHOICES, SS_CHOICES
from pgd_search.models import saveSearchForm
from pgd_splicer.chi import CHI_MAP, PROTEIN_ORDER
from pgd_splicer.sidechain import *
from django.core import serializers

#This might be a big faux pas
from pgd_splicer.models import pdb_select_settings

json_sidechain_lengths_lookup = json.dumps(bond_lengths_string_dict)
json_sidechain_angles_lookup = json.dumps(bond_angles_string_dict)



def search(request):
    """
    Handler for search form.
    """
    if request.method == 'POST': # If the form has been submitted
        form = SearchForm(request.POST) # A form bound to the POST data
        if form.is_valid(): # All validation rules pass
            #process search form into search object, remove any properties
            #that do not have values.
            data = form.cleaned_data
            for key in filter(lambda x: data[x]==None or data[x] == '', data):
                del data[key]
            search_object = Search()
            search_object.data = data
            
            #at least for now limit the size of the result set
            count = search_object.querySet().count()
            if count > settings.QUERY_LIMIT:
                form._errors['Result Size'] = ErrorList(['Your query returned more than 20,000 records, refine your search'])
            
            elif count == 0 and False:
                form._errors['Result Size'] = ErrorList(['Your query returned no results'])
            
            else:
                #store search in session
                
                search.dataset_version = pdb_select_settings.DATA_VERSION
                request.session['search'] = pickle.dumps(search_object)
                return redirect('%s/search/results/' % settings.SITE_ROOT) # Redirect after POST
        
        # package aa_choices and ss_choices
        aa_choices = []
        ss_choices = []
        for i in RESIDUE_INDEXES:
            aa_chosen = request.POST.getlist('aa_%i' % i)
            aa_choices.append([(c[0],c[1],'checked' if c[0] in aa_chosen else '') for c in AA_CHOICES])
            ss_chosen = request.POST.getlist('ss_%i' % i)
            ss_choices.append([(c[0],c[1],'checked' if c[0] in ss_chosen else '') for c in SS_CHOICES])
    else:
        aa_choices = [AA_CHOICES for i in range(settings.SEGMENT_SIZE)]
        ss_choices = [SS_CHOICES for i in range(settings.SEGMENT_SIZE)]
        form = SearchForm() # An unbound form

    #delete session variable
    try:
        del request.session['search']
    except:
        pass
    #construct a list of values for i
    iValues = [
                    # Generates a series of tuples (<value>,<signed string of value>); zero is 'i'
                    (i,'%+i'%i if i else 'i') for i in RESIDUE_INDEXES
              ]

    #order the residue properties in way that django template can handle it better
    residueFields = []
    fields = ["ss", "aa", "phi", "psi", "ome", "omep", "chi1", "chi2", "chi3", "chi4", "chi5", "bm", "bs", "bg", "h_bond_energy", "zeta", 'a1','a2','a3','a4','a5','a6','a7','L1','L2','L3','L4','L5']
    fields += sidechain_length_relationship_list
    fields += sidechain_angle_relationship_list
    for i in RESIDUE_INDEXES:
        dict = {}
        for prefix in fields:
            dict[prefix] =  form['%s_%i' % (prefix, i)]
            dict['%s_i' % prefix] =  form['%s_i_%i' % (prefix, i)]
            dict['index'] = i
        residueFields.append(dict)

    return render_to_response('search.html', {
        'form': form,
        'maxLength' : settings.SEGMENT_SIZE,
        'iValues':iValues,
        'residueFields':residueFields,
        'aa_choices':aa_choices,
        'ss_choices':ss_choices,
        'sidechain_angle_list':sidechain_angle_relationship_list,
        'sidechain_length_list':sidechain_length_relationship_list,
        'sidechain_length_lookup':json_sidechain_lengths_lookup,
        'sidechain_angle_lookup':json_sidechain_angles_lookup
    }, context_instance=RequestContext(request, processors=[settings_processor]))


def residue_is_selected(data, index):
    """
    Determine whether a given residue index is selected in search data.

    This is probably useful enough to live somewhere else, but I just don't
    fucking care right now. Sorry.
    """

    count = int(data["residues"])

    bottom = 0 - (count - 1) / 2
    top  = int(math.ceil((count - 1) / 2.0))+1
    indices = range(bottom, top)

    return index in indices


def editSearch(request, search_id=None):
    """
    Handler for editing an existing search
    """

    #load the search passed in
    if search_id:
        search = Search.objects.get(id=search_id)
        if search.user != request.user and search.isPublic == False:
            return HttpResponse("<p style='text-align:center;'>You don't have access to this search</p>")

    #else use the search in the session if it exists
    else:
        try:
            search = pickle.loads(request.session['search'])
        except KeyError:
            search = None

    if search:
        # Hax. Prepare a dataset which contains the initial angles for any
        # residues which might not have been selected in the previous search, and
        # then put the previous search's data over that dataset. This effectively
        # is the same as creating the form unbound with initial data, and then
        # binding it to the new dataset, but Django forms aren't capable of doing
        # this.
        # This fixes #1565 and related things, and could go away if the forms are
        # refactored to use FormSets and so forth.
        data = {}
        for i in RESIDUE_INDEXES:
            if not residue_is_selected(search.data, i):
                data["ome_%d" % i] = "<=-90,>=90"
                data["bm_%d" % i] = "<25"
                data["bg_%d" % i] = "<25"
                data["bs_%d" % i] = "<25"
        data.update(search.data)
        form = SearchForm(data)
    else:
        form = SearchForm() # An unbound form


    #construct a list of values for i
    iValues = [
                    # Generates a series of tuples (<value>,<signed string of value>); zero is 'i'
                    (i,'%+i'%i if i else 'i') for i in RESIDUE_INDEXES
              ]

    #order the residue properties in way that django template can handle it better
    residueFields = []
    aa_choices = []
    ss_choices = []
    fields = ["ss", "aa", "phi", "psi", "ome", "chi1", "chi2", "chi3", "chi4", "chi5", "bm", "bs", "bg", "h_bond_energy", "zeta", 'a1','a2','a3','a4','a5','a6','a7','L1','L2','L3','L4','L5']
    fields += sidechain_length_relationship_list
    fields += sidechain_angle_relationship_list
    for i in RESIDUE_INDEXES:
        valid = form.is_valid()
        aa_chosen = form.cleaned_data['aa_%i' % i]
        aa_choices.append([(c[0],c[1],'checked' if c[0] in aa_chosen else '') for c in AA_CHOICES])
        ss_chosen = form.cleaned_data['ss_%i' % i]
        ss_choices.append([(c[0],c[1],'checked' if c[0] in ss_chosen else '') for c in SS_CHOICES])

        dict = {}
        for prefix in fields:
            dict[prefix] =  form['%s_%i' % (prefix, i)]
            dict['%s_i' % prefix] =  form['%s_i_%i' % (prefix, i)]
            dict['index'] = i
        residueFields.append(dict)
    
    
    return render_to_response('search.html', {
        'form': form,
        'maxLength' : settings.SEGMENT_SIZE,
        'iValues':iValues,
        'residueFields':residueFields,
        'aa_choices':aa_choices,
        'ss_choices':ss_choices,
        'sidechain_angle_list':sidechain_angle_relationship_list,
        'sidechain_length_list':sidechain_length_relationship_list,
        'sidechain_length_lookup':json_sidechain_lengths_lookup,
        'sidechain_angle_lookup':json_sidechain_angles_lookup
    }, context_instance=RequestContext(request, processors=[settings_processor]))


def protein_search(request):
    """
    handler for ajax search box that looks up proteins.  Search box
    expects a query that is a comma delimited list of proteins.
    The last protein is autocompleted
    """
    limit = request.GET['limit']
    query = request.GET['q']

    # parse out just the last protein
    code = query.split(',')[-1].strip()

    results = Protein.objects.filter(code__istartswith=code).values_list('code')
    if len(results):
        # there were results, build a newline delimeted response as required
        # per jquery autocomplete script
        response = ''
        for issue in results:
            response = '%s\n%s' % (response, issue[0])

        return HttpResponse(response)
    else:
        return HttpResponse('')


def saved(request):

    searches = Search.objects.filter(user=request.user)

    paginator = Paginator(searches, 200000) # Show enough searches per page to effectively disable this feature for now

    # Make sure page request is an int. If not, deliver first page.
    try:
        page = int(request.GET.get('page', '1'))
    except ValueError:
        page = 1

    # If page request (9999) is out of range, deliver last page of results.
    try:
        paginatedSearch = paginator.page(page)
    except (EmptyPage, InvalidPage):
        paginatedSearch = paginator.page(paginator.num_pages)

    return render_to_response('saved.html', {
        'searches': paginatedSearch,
    }, context_instance=RequestContext(request,processors=[settings_processor]))

def help(request):
    return render_to_response('help.html', context_instance=RequestContext(request,processors=[settings_processor]))

def chi_help(request):
    return render_to_response('chi_help.html', context_instance=RequestContext(request,CHI_MAP,PROTEIN_ORDER,processors=[settings_processor]))

def qtiphelp(request):
    return render_to_response('qtiphelp.html')
    
def saveSearch(request,search_id=None):
    if request.method == 'POST':
        form = saveSearchForm(request.POST)
        if form.is_valid():
            if search_id:
                search = Search.objects.get(id=search_id)
                if request.user!=search.user:
                    return HttpResponse("<p style='text-align:center;'>You don't have access to this search</p>")
            else:
                search = pickle.loads(request.session['search'])

            data = form.cleaned_data
            search.title = data['title']
            search.description = data['description']
            search.user=request.user
            search.timestamp=datetime.now()
            search.isPublic = data['isPublic']
            search.save()

            return HttpResponseRedirect('%s/search/saved/' % settings.SITE_ROOT)
    else:
        if search_id:
            oldsearch = Search.objects.get(id=search_id)
            if request.user!=oldsearch.user:
                return HttpResponse("<p style='text-align:center;'>You don't have access to this search</p>")
            form = saveSearchForm({'title':oldsearch.title,'description':oldsearch.description,'isPublic':oldsearch.isPublic, 'search_id':search_id})
        else:
            form = saveSearchForm()
    return render_to_response('saveSearch.html', {'form': form },context_instance=RequestContext(request,processors=[settings_processor]))

def deleteSearch(request,search_id=None):
    if search_id:
        search = Search.objects.get(id=search_id)
        if search.user != request.user:
            return HttpResponse("<p style='text-align:center;'>You don't have access to this search</p>")
        search.delete()
        return HttpResponseRedirect('%s/search/saved' % settings.SITE_ROOT)
    else:
        return HttpResponseRedirect('%s/search/saved' % settings.SITE_ROOT)
