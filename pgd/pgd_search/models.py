from django.db import models
from django.contrib.auth.models import User
from pgd.pgd_core.models import Protein,Residue
from pgd.constants import AA_CHOICES, SS_CHOICES, SEQUENCE_SIZE, Subscripter

# Search
# A search query submitted by a user
class Search(models.Model):
	user = models.ForeignKey(User)

# Search_code
# Codes for the proteins searched on
class Search_code(models.Model):
	search = models.ForeignKey(Search, related_name='codes')
	code = models.CharField(max_length=4)

# Search_residue
# The search fields per residue
class Search_residue(models.Model):
	search = models.ForeignKey(Search, related_name='residues')
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
	bm = models.CharField(max_length=30)
	bs = models.CharField(max_length=30)
	bg = models.CharField(max_length=30)
	h_bond_energy = models.CharField(max_length=30)
	zeta = models.CharField(max_length=30)
	terminal_flag = models.BooleanField()
	xpr = models.BooleanField() # this field may not be necessary; it has never been implemented

# A base class for the Sequence object
class Sequence_abstract(models.Model):
	
	newID = models.ForeignKey(Protein)
	chainID = models.CharField(max_length=1)
       
	def __init__(self):
		models.Model.__init__(self)

		# make the Residue objects subscriptable
		Subscripter('residue_object', self)

	# make the residue data subscriptable? (todo?)
#	def __getattr__(self, key):
#		if key not in self.__dict__.keys() and key.startswith("residue_data_"):
#			
#		object.__getattr__(self,key)
	
	class Meta:
		abstract = True

# Build a dict for the fields of variable number
seq_dict = {'__module__' : 'pgd.pgd_search.models'}
for i in range(SEQUENCE_SIZE):
	
	# Allow access to the master Residue object...
	seq_dict["residue_object_%i" % i] = models.ForeignKey(Residue)

	# ...but make its data available (for filters, etc.) without instantiation
	seq_dict["r%i_index" % i] = models.PositiveIntegerField()
	seq_dict["r%i_newID" % i] = models.IntegerField()
	seq_dict["r%i_oldID" % i] = models.CharField(max_length=5)
	seq_dict["r%i_ss"] = models.CharField(max_length=1, choices=SS_CHOICES)
	seq_dict["r%i_terminal_flag" % i] = models.BooleanField()
	seq_dict["r%i_xpr" % i] = models.BooleanField() # probably should be replaced

	# the loops here are just to save on typing
	for j in range(1,8):
		seq_dict["r%i_a%i" % (i,j)] = models.FloatField()
	for j in range(1,6):
		seq_dict["r%i_L%i" % (i,j)] = models.FloatField()
	for j in ("phi", "psi", "ome", "chi", "bm", "bs", "bg", "h_bond_energy", "zeta"):
		seq_dict["r%i_%s" % (i,j)] = models.FloatField()

# Create the Sequence object with the fields from the dict
Sequence = type('Sequence', (Sequence_abstract,), seq_dict)
