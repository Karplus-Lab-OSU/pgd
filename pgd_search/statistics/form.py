from django import forms

from pgd_search.views import RESIDUE_INDEXES

class StatsForm(forms.Form):
    """
    Form used by the plotting function
    """
    index = forms.ChoiceField(choices=[(i,'i') if i == 0 else (i,i) for i in RESIDUE_INDEXES], initial=0)
