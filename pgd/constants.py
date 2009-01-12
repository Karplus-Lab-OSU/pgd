# The possible values for the 'aa' field (in Protein and elsewhere)
AA_CHOICES = (
	('a', 'ALA'),
	('r', 'ARG'),
	('n', 'ASN'),
	('d', 'ASP'),
    ('c', 'CYS'),
	('q', 'GLN'),
	('e', 'GLU'),
	('g', 'GLY'),
	('h', 'HIS'),
	('i', 'ILE'),
	('l', 'LEU'),
	('k', 'LYS'),
	('m', 'MET'),
	('f', 'PHE'),
	('p', 'PRO'),
	('s', 'SER'),
	('t', 'THR'),
	('w', 'TRP'),
	('y', 'TYR'),
	('v', 'VAL'),
)

# The possible values for the 'ss' field (in Protein and elsewhere)
SS_CHOICES = (
	('G', '3-turn helix (internal)'),
	('g', '3-turn helix (terminal)'),
	('H', '4-turn alpha helix (internal)'),
	('h', '4-turn alpha helix (terminal)'),
	('I', '5-turn pi helix (internal)'),
	('i', '5-turn pi helix (terminal)'),
	('T', 'H-bonded turn (internal)'),
	('t', 'H-bonded turn (terminal)'),
	('E', 'Beta sheet (internal)'),
	('e', 'Beta sheet (terminal)'),
	('B', 'Isolated beta-bridge'),
	('S', 'Bend'),
)

PLOT_PROPERTY_CHOICES = [
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
