from django.db import models

# Choice
AA_CHOICES = (
	('g' = 'gly'),
	('a' = 'ala'),
	('v' = 'val'),
	('l' = 'leu'),
	('i' = 'ile'),
	('m' = 'met'),
	('f' = 'phe'),
	('w' = 'trp'),
	('p' = 'pro'),
	('s' = 'ser'),
	('t' = 'thr'),
	('c' = 'cys'),
	('y' = 'tyr'),
	('n' = 'asn'),
	('q' = 'gln'),
	('d' = 'aps'),
	('e' = 'glu'),
	('k' = 'lys'),
	('r' = 'arg'),
	('h' = 'his'),
)
SS_TYPES = (
	('G' = '3-turn helix (internal)'),
	('g' = '3-turn helix (terminal)'),
	('H' = '4-turn alpha helix (internal)'),
	('h' = '4-turn alpha helix (terminal)'),
	('I' = '5-turn pi helix (internal)'),
	('i' = '5-turn pi helix (terminal)'),
	('T' = 'H-bonded turn (internal)'),
	('t' = 'H-bonded turn (terminal)'),
	('E' = 'Beta sheet (internal)'),
	('e' = 'Beta sheet (terminal)'),
	('B' = 'Isolated beta-bridge'),
	('S' = 'bend'),
)

# residue model
#
# Note: fields that have been changed significantly are indicated (e.g. 'aa' is
# 	now stored differently to conserve space). Also the 'idx' field is not
# 	present here because it was never implemented in the old table.
class residue(models.Model):
	code = models.CharField(max_length=4)
	aa = models.CharField(max_length=1, choices=AA_CHOICES) # new type
	chainID = models.CharField(max_length=1) 
	id = models.IntegerField()
	oldID = models.CharField(max_length=5) # is there something weird going on here?
	a1 = models.FloatField()
	a2 = models.FloatField()
	a3 = models.FloatField()
	a4 = models.FloatField()
	a5 = models.FloatField()
	a6 = models.FloatField()
	a7 = models.FloatField()
	L1 = models.FloatField()
	L2 = models.FloatField()
	L3 = models.FloatField()
	L4 = models.FloatField()
	L5 = models.FloatField()
	ss = models.CharField(max_length=1, choice=SS_TYPES) # new type (was blob, but all entries 1 char)
	phi = models.FloatField()
	psi = models.FloatField()
	ome = models.FloatField()
	chi = models.FloatField()
	ome = models.FloatField()
	bm = models.FloatField()
	bs = models.FloatField()
	H_bond_energy = models.FloatField()
	Zeta = models.FloatField()
	terminal_flag = models.BooleanField() 
	xpr = models.BooleanField() # this field may not be necessary; it has never been implemented

class residue(models.Model):

MAXNUM=10

class protein(models.Model): # incomplete
	code = models.CharField(max_length=4)
	threshold = models.IntegerField()
	resolution = models.FloatField()
	rfactor = models.FloatField()
)


#for ():
	#variable_i = djangodatatype
