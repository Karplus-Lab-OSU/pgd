from django.http import HttpResponse, HttpResponseRedirect
from django.template import RequestContext
from django.shortcuts import render_to_response
from django.conf import settings
from django.forms.util import ErrorList
from django.core.paginator import Paginator, InvalidPage, EmptyPage

import math

from pgd_core.models import Protein
from pgd_search.models import Search, Search_residue, Search_code, searchSettings
from pgd_search.views import RESIDUE_INDEXES, settings_processor
from SearchForm import SearchSyntaxField, SearchForm
from pgd_constants import AA_CHOICES, SS_CHOICES

import re

"""
Handler for search form.
"""
def search(request):

    if request.method == 'POST': # If the form has been submitted
        form = SearchForm(request.POST) # A form bound to the POST data
        if form.is_valid(): # All validation rules pass
            #process search form into search object
            search = processSearchForm(form)

            #at least for now limit the size of the result set
            if search.querySet().count() > searchSettings.query_limit:
                form._errors['Result Size'] = ErrorList(['Your query returned more than 20,000 records, refine your search'])

            elif search.querySet().count() == 0:
                form._errors['Result Size'] = ErrorList(['Your query returned no results'])

            else:
                #store search in session
                request.session['search'] = search

                return HttpResponseRedirect('%ssearch/results/' % settings.SITE_ROOT) # Redirect after POST


        # package aa_choices and ss_choices
        aa_choices = []
        ss_choices = []
        for i in RESIDUE_INDEXES:
            aa_chosen = request.POST.getlist('aa_%i' % i)
            aa_choices.append([(c[0],c[1],'checked' if c[0] in aa_chosen else '') for c in AA_CHOICES])
            ss_chosen = request.POST.getlist('ss_%i' % i)
            ss_choices.append([(c[0],c[1],'checked' if c[0] in ss_chosen else '') for c in SS_CHOICES])
    else:
        aa_choices = [AA_CHOICES for i in range(searchSettings.segmentSize)]
        ss_choices = [SS_CHOICES for i in range(searchSettings.segmentSize)]
        form = SearchForm() # An unbound form

    #construct a list of values for i
    iValues = [
                    # Generates a series of tuples (<value>,<signed string of value>); zero is 'i'
                    (i,'%+i'%i if i else 'i') for i in RESIDUE_INDEXES
              ]

    #order the residue properties in way that django template can handle it better
    residueFields = []
    for i in RESIDUE_INDEXES:
        dict = {}
        for prefix in ("ss", "aa", "phi", "psi", "ome", "chi", "bm", "bs", "bg", "h_bond_energy", "zeta", 'a1','a2','a3','a4','a5','a6','a7','L1','L2','L3','L4','L5'):
            dict[prefix] =  form['%s_%i' % (prefix, i)]
            dict['%s_i' % prefix] =  form['%s_i_%i' % (prefix, i)]
            dict['index'] = i
        residueFields.append(dict)

    return render_to_response('search.html', {
        'form': form,
        'maxLength' : searchSettings.segmentSize,
        'iValues':iValues,
        'residueFields':residueFields,
        'aa_choices':aa_choices,
        'ss_choices':ss_choices
    }, context_instance=RequestContext(request, processors=[settings_processor]))


"""
Handler for editing an existing search
"""
def editSearch(request, search_id=None):
    #load the search passed in
    if search_id:
        search = Search.objects.get(id=search_id)
        form = processSearchObject(search)


    #else use the search in the session if it exists
    else:
        try:
            search = request.session['search']
            form = processSearchObject(search)
        except KeyError:
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
    for i in RESIDUE_INDEXES:
        valid = form.is_valid()
        aa_chosen = form.cleaned_data['aa_%i' % i]
        aa_choices.append([(c[0],c[1],'checked' if c[0] in aa_chosen else '') for c in AA_CHOICES])
        ss_chosen = form.cleaned_data['ss_%i' % i]
        ss_choices.append([(c[0],c[1],'checked' if c[0] in ss_chosen else '') for c in SS_CHOICES]) 

        dict = {}
        for prefix in ("aa", "ss", "phi", "psi", "ome", "chi", "bm", "bs", "bg", "h_bond_energy", "zeta", 'a1','a2','a3','a4','a5','a6','a7','L1','L2','L3','L4','L5'):
            dict[prefix] =  form['%s_%i' % (prefix, i)]
            dict['%s_i' % prefix] =  form['%s_i_%i' % (prefix, i)]
            dict['index'] = i
        residueFields.append(dict)

    return render_to_response('search.html', {
        'form': form,
        'maxLength' : searchSettings.segmentSize,
        'iValues':iValues,
        'residueFields':residueFields,
        'aa_choices':aa_choices,
        'ss_choices':ss_choices
    }, context_instance=RequestContext(request, processors=[settings_processor]))


"""
Encode a list of AA choices into an integer value
"""
def encodeAA(list):
    aa = 0
    for i in range(len(AA_CHOICES)):
        if AA_CHOICES[i][0] in list:
            # bitwise or to add value
            aa = aa | (1 << i)
    return aa


"""
Decode an integer into a list of AA choices
"""
def decodeAA(val):
    list = []
    for i in range(len(AA_CHOICES)):
        # bitwise shift check value
        if (val & (1 << i)) != 0:
            list.append(AA_CHOICES[i][0])

    return list

"""
Encode a list of AA choices into an integer value
"""
def encodeSS(list):
    aa = 0
    for i in range(len(SS_CHOICES)):
        if SS_CHOICES[i][0] in list:
            # bitwise or to add value
            aa = aa | (1 << i)
    return aa


"""
Decode an integer into a list of SS choices
"""
def decodeSS(val):
    list = []
    for i in range(len(SS_CHOICES)):
        # bitwise shift check value
        if (val & (1 << i)) != 0:
            list.append(SS_CHOICES[i][0])

    return list



def processSearchForm(form):
    """
    Process a search form copying its data into a search object
    """
    data = form.cleaned_data
    length = searchSettings.segmentSize

    #create a new search object
    search = Search()

    #get protein properties
    search.segmentLength = int(data['residues'])
    search.resolution_min = float(data['resolutionMin']) if data['resolutionMin'] else None
    search.resolution_max = float(data['resolutionMax']) if data['resolutionMax'] else None
    search.threshold      = int(data['threshold']) if data['threshold'] else None
    search.rfactor_min    = float(data['rfactorMin']) if data['rfactorMin'] else None
    search.rfactor_max    = float(data['rfactorMax']) if data['rfactorMax'] else None
    search.rfree_min      = float(data['rfreeMin']) if data['rfreeMin'] else None
    search.rfree_max      = float(data['rfreeMax']) if data['rfreeMax'] else None

    #save search object so its residue parameters can be added
    search.save()

    #get list of proteins to filter
    if data['proteins_i'] != None:
        search.codes_include = data['proteins_i']
        for value in filter(lambda x:x!='', data['proteins'].split(',')) :
            searchCode = Search_code()
            searchCode.code = value.strip().upper()
            search.codes.add(searchCode)

    #process per residue properties
    start = 0 - (search.segmentLength-1) / 2
    stop  = int(math.ceil((search.segmentLength-1) / 2.0))+1
    for i in range(start, stop, 1):
        hasField = False
        residue = Search_residue()
        residue.index   = i

        #process ss
        if data['ss_%i' % i]:
            residue.ss_int          = encodeSS(data['ss_%i' % i])
            residue.ss_int_include  = data['ss_i_%i' % i]
            hasField = True

        #process aa
        if data['aa_%i' % i]:
            residue.aa_int          = encodeAA(data['aa_%i' % i])
            residue.aa_int_include  = data['aa_i_%i' % i]
            hasField = True

        #process all other fields
        for prefix in ("phi", "psi", "ome", "chi", "bm", "bs", "bg", "h_bond_energy", "zeta", 'a1','a2','a3','a4','a5','a6','a7','L1','L2','L3','L4','L5'):
            if data['%s_%i' % (prefix, i)]:
                residue.__dict__[prefix] = data['%s_%i' % (prefix, i)]
                residue.__dict__['%s_include' % prefix] = data['%s_i_%i' % (prefix, i)]
                hasField = True

        # only add the residue if there was a value in the field
        if hasField:
            search.residues.add(residue)
            residue.search = search
            residue.save()

    return search


def processSearchObject(search):
    """
    Process a search object copying its values into a searchForm
    """

    data = {
        #get protein properties
        'residues'      :search.segmentLength,
        'resolutionMin' :search.resolution_min,
        'resolutionMax' :search.resolution_max,
        'rfreeMin'      :search.rfree_min,
        'rfreeMax'      :search.rfree_max,
        'rfactorMin'    :search.rfactor_min,
        'rfactorMax'    :search.rfactor_max,
        'threshold'     :search.threshold
    }

    #get list of proteins to filter
    if search.codes_include != None:
        data['proteins_i'] = search.codes_include
        codes = []
        for code in search.codes.all():
            codes.append(code.code)
        data['proteins'] = ', '.join(codes)

    #setup defaults - initial values are not set when passing a dict to the form constructor
    # so any value with a default value must be initialized prior to values are pulled out
    # of the Search object
    for i in RESIDUE_INDEXES:
        data['bm_%i' % i ]  = SearchForm.base_fields['bm_%i' % i].initial
        data['bg_%i' % i ]  = SearchForm.base_fields['bg_%i' % i].initial
        data['bs_%i' % i ]  = SearchForm.base_fields['bs_%i' % i].initial
        data['ome_%i' % i ] = SearchForm.base_fields['ome_%i' % i].initial

    #process per residue properties
    for residue in search.residues.all():
        i = residue.index

        #process ss
        if residue.ss_int:
            data['ss_%i' % i]   = decodeSS(residue.ss_int)
            data['ss_i_%i' % i] = residue.ss_int_include

        #process aa
        if residue.aa_int:
            data['aa_%i' % i]   = decodeAA(residue.aa_int)
            data['aa_i_%i' % i] = residue.aa_int_include

        #process all other fields
        for prefix in ("phi", "psi", "ome", "chi", "bm", "bs", "bg", "h_bond_energy", "zeta", 'a1','a2','a3','a4','a5','a6','a7','L1','L2','L3','L4','L5'):
            data['%s_%i' % (prefix, i)] = residue.__dict__[prefix]
            data['%s_i_%i' % (prefix, i)] = residue.__dict__['%s_include' % prefix]

    return SearchForm(data)


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

    searches = Search.objects.all()

    paginator = Paginator(searches, 20) # Show 20 searches per page

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
    }, context_instance=RequestContext(request))

def help(request):
    return render_to_response('help.html', context_instance=RequestContext(request))

def qtiphelp(request):
    return render_to_response('qtiphelp.html', context_instance=RequestContext(request))

