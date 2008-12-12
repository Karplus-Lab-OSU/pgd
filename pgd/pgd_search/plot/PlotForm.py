from django import forms

from pgd_search.views import RESIDUE_INDEX_START, RESIDUE_INDEX_STOP

#choice for occurence of property
ATTRIBUTE_CHOICES = [
                    ("Observations","Observations"),
                    ("L1","L1"),
                    ("L2","L2"),
                    ("L3","L3"),
                    ("L4","L4"),
                    ("L5","L5"),
                    ("a1","a1"),
                    ("a2","a2"),
                    ("a3","a3"),
                    ("a4","a4"),
                    ("a5","a5"),
                    ("a6","a6"),
                    ("a7","a7"),
                    ("ome","ome"),
                    ("chi","chi"),
                    ]

# choices for properties mapped to axis
PROPERTY_CHOICES = [
                    ("L1","L1"),
                    ("L2","L2"),
                    ("L3","L3"),
                    ("L4","L4"),
                    ("L5","L5"),
                    ("a1","a1"),
                    ("a2","a2"),
                    ("a3","a3"),
                    ("a4","a4"),
                    ("a5","a5"),
                    ("a6","a6"),
                    ("a7","a7"),
                    ("ome","ome"),
                    ("chi","chi"),
                    ("phi","phi"),
                    ("psi","psi"),
                    ]

"""
Form used by the plotting function
"""
class PlotForm(forms.Form):
    attribute       = forms.ChoiceField(choices=ATTRIBUTE_CHOICES, initial='Observations')
    xProperty       = forms.ChoiceField(choices=PROPERTY_CHOICES, initial='phi')
    yProperty       = forms.ChoiceField(choices=PROPERTY_CHOICES, initial='psi')
    reference       = forms.FloatField(required=False, widget=forms.TextInput(attrs={'size':8}))
    x               = forms.FloatField(initial=-180, widget=forms.TextInput(attrs={'size':4}))
    x1              = forms.FloatField(initial=180, widget=forms.TextInput(attrs={'size':4}))
    y               = forms.FloatField(initial=-180, widget=forms.TextInput(attrs={'size':4}))
    y1              = forms.FloatField(initial=180, widget=forms.TextInput(attrs={'size':4}))
    residue         = forms.ChoiceField(choices=[(i,'i') if i == 0 else (i,i) for i in range(RESIDUE_INDEX_START,RESIDUE_INDEX_STOP)], initial=0)
    xBin            = forms.FloatField(initial=10, widget=forms.TextInput(attrs={'size':4}))
    yBin            = forms.FloatField(initial=10, widget=forms.TextInput(attrs={'size':4}))