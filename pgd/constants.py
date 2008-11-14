SEQUENCE_SIZE = 10

# Choice
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

def pack_dict(label,object,limit):
	mydict = dict(("%s_%d" % (label,key), object()) for key in range(limit))
	mydict['__module__'] = 'pgd.tasks.models'
	return mydict

# subscripter is a way to make a property subscriptable
# this allows you to reference properties by index
class Subscripter():
    def __init__(self, key, parent):
        self.key = key
        self.parent = parent
        #add this instance to the parent. doing this here
        #makes defining subscriptor instance simpler because
        #you only need to specify the key once
        parent.__dict__[key] = self

    def __getitem__(self, i):
        return self.parent.__dict__['%s_%d' % (self.key, i)]

    def __setitem__(self, i, val):
        self.parent.__dict__['%s_%d' % (self.key, i)] = val
