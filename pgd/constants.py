SEQUENCE_SIZE = 10

# The possible values for the 'aa' field (in Protein and elsewhere)
AA_CHOICES = (
	('g', 'gly'),
	('a', 'ala'),
	('v', 'val'),
	('l', 'leu'),
	('i', 'ile'),
	('m', 'met'),
	('f', 'phe'),
	('w', 'trp'),
	('p', 'pro'),
	('s', 'ser'),
	('t', 'thr'),
	('c', 'cys'),
	('y', 'tyr'),
	('n', 'asn'),
	('q', 'gln'),
	('d', 'aps'),
	('e', 'glu'),
	('k', 'lys'),
	('r', 'arg'),
	('h', 'his'),
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
	('S', 'bend'),
)

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
