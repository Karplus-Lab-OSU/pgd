from django.db import models
from django.contrib.auth.models import User
from pgd.pgd_core.models import Protein
from pgd.constants import AA_CHOICES, SS_CHOICES, pack_dict, Subscripter

# residue model
# (was 'protein')
#
# Note: The 'idx' field is not present here because it was never implemented in 
# 	the old table, and the model uses a default primary key field (can't
# 	spread a key across fields in django... yet?).
#class Residue(models.Model):
#	code = models.ForeignKey(protein)
#	aa = models.CharField(max_length=1, choices=AA_CHOICES) # new type
#	chainID = models.CharField(max_length=1) 
#	id = models.IntegerField()
#	oldID = models.CharField(max_length=5) # is there something weird going on here?
#	a1 = models.FloatField()
#	a2 = models.FloatField()
#	a3 = models.FloatField()
#	a4 = models.FloatField()
#	a5 = models.FloatField()
#	a6 = models.FloatField()
#	a7 = models.FloatField()
#	L1 = models.FloatField()
#	L2 = models.FloatField()
#	L3 = models.FloatField()
#	L4 = models.FloatField()
#	L5 = models.FloatField()
#	ss = models.CharField(max_length=1, choice=SS_TYPES) # new type (was blob, but all entries 1 char)
#	phi = models.FloatField()
#	psi = models.FloatField()
#	ome = models.FloatField()
#	chi = models.FloatField()
#	ome = models.FloatField()
#	bm = models.FloatField()
#	bs = models.FloatField()
#	h_bond_energy = models.FloatField()
#	zeta = models.FloatField()
#	terminal_flag = models.BooleanField() 
#	xpr = models.BooleanField() # this field may not be necessary; it has never been implemented

# protein model
# (was 'protein_info')
#class Protein(models.Model):
#	code = models.CharField(max_length=4, primary_key=True, unique=True)
#	threshold = models.IntegerField() # new type; this should probably be a boolean type
#	resolution = models.FloatField()
#	rfactor = models.FloatField()

# sequence model
# 
# Note: uses django-magic!
#seq_dict = {
#	'code' : models.ForeignKey(protein),
#	'chainID' : models.CharField(max_length=1),
#	'id' : models.IntegerField(),
#	'oldID' : models.CharField(max_length=5), # is there something weird going on in this field?
#}
#for i in ('phi','psi','ome','chi','bm','bs','bg','H_bond_energy','zeta'):
#	seq_dict[i] = models.FloatField()
#for i in range(0,SEQUENCE_SIZE):
#	seq_dict['%d_ss'] : models.CharField(max_length=1, choice=SS_TYPES), # new type
#	seq_dict['%d_aa' % i] : models.CharField(max_length=1, choices=AA_CHOICES), # new type
#	for j in range(1,8):
#		seq_dict['%d_a%d' % (i,j)] = models.FloatField()
#	for j in range(1,6):
#		seq_dict['%d_L%d' % (i,j)] = models.FloatField()


class Search(models.Model):
	user = models.ForeignKey(User)
	code = models.ForeignKey(Protein)

class Search_code(models.Model):
	search = models.ForeignKey(Search, related_name='codes')
	code = models.CharField(max_length=4)

class Search_residue(models.Model):
	search = models.ManyToManyField(Search, related_name='residues')
	index = models.PositiveIntegerField()
	chainID = models.CharField(max_length=1)
	newID = models.IntegerField()
	oldID = models.CharField(max_length=5) # perhaps we don't need to search on this field
	a1 = models.CharField(max_length=30)
	a2 = models.CharField(max_length=30)
	a3 = models.CharField(max_length=30)
	a4 = models.CharField(max_length=30)
	a5 = models.CharField(max_length=30)
	a6 = models.CharField(max_length=30)
	a7 = models.CharField(max_length=30)
	L1 = models.CharField(max_length=30)
	L2 = models.CharField(max_length=30)
	L3 = models.CharField(max_length=30)
	L4 = models.CharField(max_length=30)
	L5 = models.CharField(max_length=30)
	ss = models.CharField(max_length=1, choices=SS_CHOICES) # new type (was blob, but all entries 1 char)
	phi = models.CharField(max_length=30)
	psi = models.CharField(max_length=30)
	ome = models.CharField(max_length=30)
	chi = models.CharField(max_length=30)
	ome = models.CharField(max_length=30)
	bm = models.CharField(max_length=30)
	bs = models.CharField(max_length=30)
	bs = models.CharField(max_length=30)
	h_bond_energy = models.CharField(max_length=30)
	zeta = models.CharField(max_length=30)
	terminal_flag = models.BooleanField() 
	xpr = models.BooleanField() # this field may not be necessary; it has never been implemented
	 
