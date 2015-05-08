from django import forms

from pgd_search.views import RESIDUE_INDEXES

#choice for occurence of property
ATTRIBUTE_CHOICES = [
                    ("Observations",'Observations'),
                    ("L1",u'C<sup>-1</sup>N'),
                    ("L2",u'NC<sup>&alpha;</sup>'),
                    ("L3",u'C<sup>&alpha;</sup>C<sup>&beta;</sup>'),
                    ("L4",u'C<sup>&alpha;</sup>C'),
                    ("L5",u'CO'),
                    ("a1",u'C<sup>-1</sup>NC<sup>&alpha;</sup>'),
                    ("a2",u'NC<sup>&alpha;</sup>C<sup>&beta;</sup>'),
                    ("a3",u'NC<sup>&alpha;</sup>C'),
                    ("a4",u'C<sup>&beta;</sup>C<sup>&alpha;</sup>C'),
                    ("a5",u'C<sup>&alpha;</sup>CO'),
                    ("a6",u'C<sup>&alpha;</sup>CN<sup>+1</sup>'),
                    ("a7",u'OCN<sup>+1</sup>'),
                    ("ome",u'&omega;'),
                    ("omep",u'&omega;P'),
                    ("chi1",u'&chi;<sup>1</sup>'),
                    ("chi2",u'&chi;<sup>2</sup>'),
                    ("chi3",u'&chi;<sup>3</sup>'),
                    ("chi4",u'&chi;<sup>4</sup>'),
                    ("chi5",u'&chi;<sup>5</sup>'),
                    ("phi",u'&Phi;'),
                    ("psi",u'&psi;'),
                    ('zeta',u'&zeta;'),
                    #('h_bond_energy','H Bond'),
                    ]

# choices for properties mapped to axis
PROPERTY_CHOICES = [
                    ("L1",u'C<sup>-1</sup>N'),
                    ("L2",u'NC<sup>&alpha;</sup>'),
                    ("L3",u'C<sup>&alpha;</sup>C<sup>&beta;</sup>'),
                    ("L4",u'C<sup>&alpha;</sup>C'),
                    ("L5",u'CO'),
                    ("a1",u'C<sup>-1</sup>NC<sup>&alpha;</sup>'),
                    ("a2",u'NC<sup>&alpha;</sup>C<sup>&beta;</sup>'),
                    ("a3",u'NC<sup>&alpha;</sup>C'),
                    ("a4",u'C<sup>&beta;</sup>C<sup>&alpha;</sup>C'),
                    ("a5",u'C<sup>&alpha;</sup>CO'),
                    ("a6",u'C<sup>&alpha;</sup>CN<sup>+1</sup>'),
                    ("a7",u'OCN<sup>+1</sup>'),
                    ("ome",u'&omega;'),
                    ("omep",u'&omega;P'),
                    ("chi1",u'&chi;<sup>1</sup>'),
                    ("chi2",u'&chi;<sup>2</sup>'),
                    ("chi3",u'&chi;<sup>3</sup>'),
                    ("chi4",u'&chi;<sup>4</sup>'),
                    ("chi5",u'&chi;<sup>5</sup>'),
                    ("phi",u'&Phi;'),
                    ("psi",u'&psi;'),
                    ('zeta',u'&zeta;'),
                    #('h_bond_energy','H Bond'),
                    ]


PROPERTY_CHOICES_DICT = {}
for prop, label in PROPERTY_CHOICES:
    PROPERTY_CHOICES_DICT[prop] = label          


BACKGROUND_CHOICES = [
                    ('#ffffff','White'),
                    ('#000000','Black'),
                    ('#666666','Gray'),
                    ('#222222','Dark Gray'),
                    (None,'Transparent'),
]

GRAPH_CHOICES = [
                    ('#222222','Dark Gray'),
                    ('#666666','Gray'),
                    ('#000000','Black'),
                    ('#ffffff','White'),
                    (None,'Transparent'),
]

TEXT_CHOICES = [
                    ('#000000','Black'),
                    ('#ffffff','White'),
                    ('#666666','Gray'),
                    ('#222222','Dark Gray'),
]

HUE_CHOICES = [
                    ('green','Green'),
                    ('blue','Blue'),
                    ('red','Red'),
                    ('black','Black/White'),
]

HASH_CHOICES = [
                    ('#666666','Gray'),
                    ('#222222','Dark Gray'),
                    ('#000000','Black'),
                    ('#ffffff','White'),
]


"""
Form used by the plotting function
"""
class PlotForm(forms.Form):
    attribute       = forms.CharField(initial='Observations')

    xProperty       = forms.CharField(initial='phi')
    yProperty       = forms.CharField(initial='psi')

    reference       = forms.FloatField(required=False,
                                       widget=forms.TextInput(attrs={'size':8}))
    sigmaVal        = forms.FloatField(initial=3,
                                       min_value=0,
                                       required=False,
                                       widget=forms.TextInput(attrs={'size':8}))

    x               = forms.FloatField(required=False,
                                       initial=-180,
                                       widget=forms.TextInput(attrs={'size':4}))
    x1              = forms.FloatField(required=False,
                                       initial=180,
                                       widget=forms.TextInput(attrs={'size':4}))

    y               = forms.FloatField(required=False,
                                       initial=-180,
                                       widget=forms.TextInput(attrs={'size':4}))
    y1              = forms.FloatField(required=False,
                                       initial=180,
                                       widget=forms.TextInput(attrs={'size':4}))

    residue_attribute = forms.ChoiceField(choices=[(i,'i') if i == 0 else (i,i) for i in RESIDUE_INDEXES],
                                          initial=0)

    residue_xproperty = forms.ChoiceField(choices=[(i,'i') if i == 0 else (i,i) for i in RESIDUE_INDEXES],
                                          initial=0)
    residue_yproperty = forms.ChoiceField(choices=[(i,'i') if i == 0 else (i,i) for i in RESIDUE_INDEXES],
                                          initial=0)

    xBin            = forms.FloatField(required=False,
                                       initial=10,
                                       widget=forms.TextInput(attrs={'size':4}))
    yBin            = forms.FloatField(required=False,
                                       initial=10,
                                       widget=forms.TextInput(attrs={'size':4}))

    #custom plot properties
    background_color= forms.ChoiceField(choices=BACKGROUND_CHOICES)
    graph_color     = forms.ChoiceField(choices=GRAPH_CHOICES)
    text_color      = forms.ChoiceField(choices=TEXT_CHOICES)
    plot_hue        = forms.ChoiceField(choices=HUE_CHOICES)
    hash_color      = forms.ChoiceField(choices=HASH_CHOICES)
    height          = forms.IntegerField(initial=470,
                                         widget=forms.TextInput(attrs={'size':4}))
    width           = forms.IntegerField(initial=560,
                                         widget=forms.TextInput(attrs={'size':4}))

    def clean(self):
        data = self.cleaned_data
        try:
            data['xProperty'] = data['xProperty'].replace('-','_')
            data['yProperty'] = data['yProperty'].replace('-','_')
            data['attribute'] = data['attribute'].replace('-','_')
        except KeyError:
            pass
        return data
