from django.db import models
from pgd.constants import AA_CHOICES, SS_CHOICES, Subscripter

# Protein model
# (was 'protein_info')
# Contains the general information about a protein
class Protein(models.Model):
	code = models.CharField(max_length=4, primary_key=True, unique=True)
	threshold = models.IntegerField() # new type; this should probably be a boolean type
	resolution = models.FloatField()
	rfactor = models.FloatField()

# Residue model
# (was 'protein')
# Contains information regarding each amino acid and its geometry
# (Note: fields need to be commented)
class Residue(models.Model):
	code = models.ForeignKey(Protein, related_name='residues')
	aa = models.CharField(max_length=1, choices=AA_CHOICES) # new type
	chainID = models.CharField(max_length=1) 
	newID = models.IntegerField()
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
	ss = models.CharField(max_length=1, choices=SS_CHOICES) # new type (was blob, but all entries 1 char)
	phi = models.FloatField()
	psi = models.FloatField()
	ome = models.FloatField()
	chi = models.FloatField()
	bm = models.FloatField()
	bs = models.FloatField()
	bg = models.FloatField()
	h_bond_energy = models.FloatField()
	zeta = models.FloatField()
	terminal_flag = models.BooleanField() 
	xpr = models.BooleanField() # this field may not be necessary; it has never been implemented
