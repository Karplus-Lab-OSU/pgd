# The possible values for the 'aa' field (in Protein and elsewhere)
AA_CHOICES = (
    ('a', 'Ala'),
    ('r', 'Arg'),
    ('n', 'Asn'),
    ('d', 'Asp'),
    ('c', 'Cys'),
    ('q', 'Gln'),
    ('e', 'Glu'),
    ('g', 'Gly'),
    ('h', 'His'),
    ('i', 'Ile'),
    ('l', 'Leu'),
    ('k', 'Lys'),
    ('m', 'Met'),
    ('f', 'Phe'),
    ('p', 'Pro'),
    ('s', 'Ser'),
    ('t', 'Thr'),
    ('w', 'Trp'),
    ('y', 'Tyr'),
    ('v', 'Val'),
)

AA_CHOICES_DICT = {}
for choice in AA_CHOICES:
    AA_CHOICES_DICT[choice[0]] = choice[1]

# The possible values for the 'ss' field (in Protein and elsewhere)
'''SS_CHOICES = (
    ('G', 'G'),
    ('g', 'g'),
    ('H', 'H'),
    ('h', 'h'),
    ('I', 'I'),
    ('i', 'i'),
    ('T', 'T'),
    ('t', 't'),
    ('E', 'E'),
    ('e', 'e'),
    ('B', 'B'),
    ('S', 'S'),
)'''
SS_CHOICES = (
    ('H', '&alpha; helix'),
    #('h', '&alpha; helix (terminal)'),
    ('G', '3<sub>10</sub> helix'),
    #('g', '3<sub>10</sub> helix (terminal)'),
    ('E', '&beta; sheet'),
    #('e', '&beta; sheet (terminal)'),
    ('T', 'Turn'),
    #('t', 'H-bonded turn (terminal)'),
    ('S', 'Bend'),
    ('B', '&beta;-bridge'),
    ('I', '&pi; helix'),
    #('i', '&pi; helix'),
)

PLOT_PROPERTY_CHOICES = [
    ("phi","phi"),
    ("psi","psi"),    
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
    ("chi1","chi1"),
    ("chi2","chi3"),
    ("chi3","chi3"),
    ("chi4","chi4"),
    ("chi5","chi5"),
    ("zeta","zeta")
    #("h_bond_energy","H Bond"),
    ]

# This class makes a property subscriptable
# (may get moved to pgd_search/models.py)
class Subscripter():
    def __init__(self, key, parent):
        self.key = key
        self.parent = parent
        #add this instance to the parent. doing this here
        #makes defining subscriptor instance simpler because
        #you only need to specify the key once
        parent.__dict__[key] = self

    def __getitem__(self, i):
        return self.parent.__dict__['%s_%i' % (self.key, i)]

    def __setitem__(self, i, val):
        self.parent.__dict__['%s_%i' % (self.key, i)] = val
