from django.http import HttpResponseRedirect
from django.template import RequestContext
from django.shortcuts import render_to_response
from django.conf import settings
from django.forms.util import ErrorList

import math

from pgd_search.models import Search, Search_residue, Search_code, searchSettings
from pgd_search.views import RESIDUE_INDEXES
from SearchForm import SearchForm
from constants import AA_CHOICES

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
            if search.querySet().count() > 10000:
                form._errors['Result Size'] = ErrorList(['Your query returned more than 10,000 records, refine your search'])

            else:
                #store search in session
                request.session['search'] = search

                return HttpResponseRedirect('%ssearch/results/' % settings.SITE_ROOT) # Redirect after POST
    else:
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
        'residueFields':residueFields
    }, context_instance=RequestContext(request))


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
Process a search form copying its data into a search object
"""
def processSearchForm(form):
    data = form.cleaned_data
    length = searchSettings.segmentSize

    #create a new search object
    search = Search()

    #get protein properties
    search.residueCount  = int(data['residues'])
    search.resolutionMin = float(data['resolutionMin'])
    search.resolutionMax = float(data['resolutionMax'])

    #save search object so its residue parameters can be added
    search.save()

    #get list of proteins to filter
    for value in data['proteins']:
        searchCode = Search_code()
        searchCode.code = value
        search.codes.add(searchCode)

    #process per residue properties
    start = 0 - (search.residueCount-1) / 2
    stop  = int(math.ceil((search.residueCount-1) / 2.0))+1
    for i in range(start, stop, 1):
        hasField = False
        residue = Search_residue()
        residue.index   = i

        #process ss
        if data['ss_%i' % i]:
            residue.ss          = data['ss_%i' % i]
            residue.ss_include  = data['ss_i_%i' % i]
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
