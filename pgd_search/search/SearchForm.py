import re

from django.conf import settings
from django import forms
from django.forms import (Field, CharField, ChoiceField, FloatField,
                          IntegerField, MultipleChoiceField, HiddenInput,
                          SelectMultiple, TextInput, Form)

from pgd_search.views import RESIDUE_INDEXES
from pgd_constants import AA_CHOICES, SS_CHOICES
from pgd_splicer.sidechain import sidechain_angle_relationship_list, sidechain_length_relationship_list

"""
Custom Field for fields that support query syntax parsing
"""
class SearchSyntaxField(Field):
    #Validates a field to make sure that it has valid syntax for a search field
    #floatNum    = r'(-?(([1-9]\d*|0)(\.\d+)?|(\.\d+)))'
    #compOp      = r'([<>]=?)'
    #argEntry    = r'(floatNum(-floatNum)?|compOpfloatNum)'
    #synPat      = r'(argEngry)(,argEntry)*'
    syntaxPattern = re.compile(
        r'^(-?(([1-9]\d*|0)(\.\d+)?|(\.\d+))(--?(([1-9]\d*|0)(\.\d+)?|(\.\d+)))?|[<>]=?-?(([1-9]\d*|0)(\.\d+)?|(\.\d+)))(,(-?(([1-9]\d*|0)(\.\d+)?|(\.\d+))(--?(([1-9]\d*|0)(\.\d+)?|(\.\d+)))?|[<>]=?-?(([1-9]\d*|0)(\.\d+)?|(\.\d+))))*$'
    )

    def clean(self, value):
        if value == None:
            return None

        if value == '':
            return value

        # remove all whitespace
        cleaned = value.replace(' ','')

        # replace or synonyms with comma
        p = re.compile( '(or|\|\|)')
        cleaned = p.sub(',', cleaned)

        # check that the field matches the syntax
        if self.syntaxPattern.match(cleaned) == None:
            raise forms.ValidationError('Input is not valid query syntax')

        return cleaned

"""
Search form used by search handler.
"""
class SearchFormBase(Form):
    threshold       = ChoiceField(choices=[(25,25),(90,90)], required=False)
    resolutionMin   = FloatField(required=False,
                                 min_value=0,
                                 initial=0,
                                 widget=TextInput(attrs={'size':3}))
    resolutionMax   = FloatField(required=False,
                                 min_value=0,
                                 initial=1.2,
                                 widget=TextInput(attrs={'size':3}))
    rfactorMin      = FloatField(required=False,
                                 min_value=0,
                                 initial=0,
                                 widget=TextInput(attrs={'size':3}))
    rfactorMax      = FloatField(required=False,
                                 min_value=0,
                                 initial=0.25,
                                 widget=TextInput(attrs={'size':3}))
    rfreeMin        = FloatField(required=False,
                                 min_value=0,
                                 initial=0,
                                 widget=TextInput(attrs={'size':3}))
    rfreeMax        = FloatField(required=False,
                                 min_value=0,
                                 initial=0.30,
                                 widget=TextInput(attrs={'size':3}))
    proteins        = CharField(required=False)
    proteins_i      = IntegerField(required=False,
                                   widget=HiddenInput(attrs={'class':'include'}))
    residues        = ChoiceField(choices=[(i,i) for i in range(1, settings.SEGMENT_SIZE+1)],
                                  initial=3)

# Build a dict for the fields of variable number
form_dict = {'__module__' : 'pgd_search.views'}

for i in RESIDUE_INDEXES:
    form_dict["aa_%i" % i]      = MultipleChoiceField(choices=AA_CHOICES,
                                                      required=False,
                                                      widget=SelectMultiple(attrs={'class':'field'}))
    form_dict["aa_i_%i" % i]    = IntegerField(required=False,
                                               widget=HiddenInput(attrs={'class':'include'}))

    form_dict["ss_%i" % i]      = MultipleChoiceField(choices=SS_CHOICES,
                                                      required=False,
                                                      widget=SelectMultiple(attrs={'class':'field'}))
    form_dict["ss_i_%i" % i]    = IntegerField(required=False,
                                               widget=HiddenInput(attrs={'class':'include'}))

    # the loops here are just to save on space/typing
    for j in range(1,8):   # Angles
        form_dict["a%i_%i" % (j,i)]     = SearchSyntaxField(required=False,
                                                            widget=TextInput(attrs={'class':'field needs_reset', 'size':8}))
        form_dict["a%i_i_%i" % (j,i)]   = IntegerField(required=False,
                                                       widget=HiddenInput(attrs={'class':'include'}))

    for j in range(1,6):   # Lengths
        form_dict["L%i_%i" % (j, i)]    = SearchSyntaxField(required=False,
                                                            widget=TextInput(attrs={'class':'field needs_reset', 'size':8}))
        form_dict["L%i_i_%i" % (j, i)]  = IntegerField(required=False,
                                                       widget=HiddenInput(attrs={'class':'include'}))

    plot_property = ['phi','psi','ome','chi1','chi2','chi3','chi4','chi5','bm',
                     'bs','bg','h_bond_energy','zeta']
    for j in plot_property:
        form_dict["%s_%i" % (j, i)]     = SearchSyntaxField(required=False,
                                                            widget=TextInput(attrs={'class':'field needs_reset', 'size':8}))
        form_dict["%s_i_%i" % (j, i)]   = IntegerField(required=False,
                                                       widget=HiddenInput(attrs={'class':'include'}))


    form_dict["ome_%i" % i]             = SearchSyntaxField(initial='<=-90,>=90',
                                                            required=False,
                                                            widget=TextInput(attrs={'class':'field needs_reset', 'size':8}))
    form_dict["ome_i_%i" % i]           = IntegerField(required=False,
                                                       widget=HiddenInput(attrs={'class':'include'}))

    form_dict["omep_%i" % i]             = SearchSyntaxField(initial='<=-90,>=90',
                                                             required=False,
                                                             widget=TextInput(attrs={'class':'field needs_reset', 'size':8}))
    form_dict["omep_i_%i" % i]           = IntegerField(required=False,
                                                        widget=HiddenInput(attrs={'class':'include'}))

    for j in ("bm", "bs", "bg"):
        form_dict["%s_%i" % (j, i)]     = SearchSyntaxField(initial='<25',
                                                            required=False,
                                                            widget=TextInput(attrs={'class':'field needs_reset', 'size':8}))
        form_dict["%s_i_%i" % (j, i)]   = IntegerField(required=False,
                                                       widget=HiddenInput(attrs={'class':'include'}))

    #Adding the occupancy into search page
    # issue : https://code.osuosl.org/issues/17577
    for j in ("occm", "occscs") :
        # form_dict["%s_%i" % (j, i)]     = FloatField(max_value=1.0, min_value=0.0, required=False)
        form_dict["%s_%i" % (j, i)]     = SearchSyntaxField(required=False,
                                                            widget=TextInput(attrs={'class':'field needs_reset', 'size':8}))
        form_dict["%s_i_%s" % (j, i)]   = IntegerField(required=False, widget=HiddenInput(attrs={'class':'include'})) 

    for j in sidechain_angle_relationship_list:
        form_dict["%s_%i" % (j,i)]      = SearchSyntaxField(required=False,
                                                            widget=TextInput(attrs={'class':'field needs_reset', 'size':8}))
        form_dict["%s_i_%i" % (j, i)]   = IntegerField(required=False,
                                                       widget=HiddenInput(attrs={'class':'include'}))

    for j in sidechain_length_relationship_list:
        form_dict["%s_%i" % (j,i)]      = SearchSyntaxField(required=False,
                                                            widget=TextInput(attrs={'class':'field needs_reset', 'size':8}))
        form_dict["%s_i_%i" % (j, i)]   = IntegerField(required=False,
                                                       widget=HiddenInput(attrs={'class':'include'}))

# Create the Search Form with the fields from the dict
SearchForm = type('SearchForm', (SearchFormBase,), form_dict)
