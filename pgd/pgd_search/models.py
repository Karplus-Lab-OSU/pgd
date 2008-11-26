from django.db import models
from django.contrib.auth.models import User
from pgd_core.models import Protein,Residue
from constants import AA_CHOICES, SS_CHOICES, SEQUENCE_SIZE, Subscripter
from exceptions import AttributeError

# Search
# A search query submitted by a user
class Search(models.Model):
    user = models.ForeignKey(User)

# Search_code
# Codes for the proteins searched on
class Search_code(models.Model):
    search  = models.ForeignKey(Search, related_name='codes')
    code    = models.CharField(max_length=4)

# Search_residue
# The search fields per residue
class Search_residue(models.Model):
    search          = models.ForeignKey(Search, related_name='residues')
    index           = models.PositiveIntegerField()
    chainID         = models.CharField(max_length=1)
    newID           = models.IntegerField()
    oldID           = models.CharField(max_length=5) # perhaps we don't need to search on this field?
    aa_int          = models.IntegerField()

    a1              = models.CharField(max_length=30)
    a2              = models.CharField(max_length=30)
    a3              = models.CharField(max_length=30)
    a4              = models.CharField(max_length=30)
    a5              = models.CharField(max_length=30)
    a6              = models.CharField(max_length=30)
    a7              = models.CharField(max_length=30)
    L1              = models.CharField(max_length=30)
    L2              = models.CharField(max_length=30)
    L3              = models.CharField(max_length=30)
    L4              = models.CharField(max_length=30)
    L5              = models.CharField(max_length=30)
    ss              = models.CharField(max_length=1, choices=SS_CHOICES) # new type (was blob, but all entries 1 char)
    phi             = models.CharField(max_length=30)
    psi             = models.CharField(max_length=30)
    ome             = models.CharField(max_length=30)
    chi             = models.CharField(max_length=30)
    bm              = models.CharField(max_length=30)
    bs              = models.CharField(max_length=30)
    bg              = models.CharField(max_length=30)
    h_bond_energy   = models.CharField(max_length=30)
    zeta            = models.CharField(max_length=30)
    terminal_flag   = models.BooleanField()
    xpr             = models.BooleanField() # this field may not be necessary; it has never been implemented

    def __init__(self):
        models.Model.__init__(self)
        # make 'residues' a subscriptable way to access Residue onjects
        Residue_subscripter(self, 'residues')
        # populate 'aa' with a dictionary of allowed 'aa' values
        self.aa = dict([(j[1],1 if self.aa_int == None else 1&self.aa_int>>i) for i,j in enumerate(AA_CHOICES)])


# Residue_subscripter
# A subscripter for iterating through the Residues in a Segment
class Residue_subscripter():
    def __init__(self, key, parent):
        self.parent = parent
        #add this instance to the parent. doing this here
        #makes defining subscriptor instance simpler because
        #you only need to specify the key once
        parent.__dict__[key] = self

    def __getitem__(self, i):
        try: # Get the object, if it's been instantiated...
            return self.parent.__dict__['r%i' % i]
        except KeyError: # ...otherwise, instantiate it.
            try: # If the object can be instantiated, return it...
                self.parent.__dict__['r%i' % i] = Residue.objects.filter(id=self.parent.__dict__['r%i_id' % i])[0]
                return self.parent.__dict__['r%i' % i]
            except KeyError: # ...otherwise, give an index error
                raise IndexError

    def __iter__(self):
        # This function makes a generator object
        def residue_generator(outer):
            for i in range(SEQUENCE_SIZE):
                try: # Get the next object...
                    yield outer.__getitem__(i)
                except IndexError: # ...until there are no more.
                    raise StopIteration
        return residue_generator(self)

# Segment_abstract
# A base class for the Segment object
class Segment_abstract(models.Model):

    code = models.ForeignKey(Protein)
    chainID = models.CharField(max_length=1)

    def __init__(self):
        models.Model.__init__(self)
        Residue_subscripter(self, 'residues')

    class Meta:
        abstract = True

# Build a dict for the fields of variable number
seq_dict = {'__module__' : 'pgd_search.models'}
for i in range(2):

    seq_dict["r%i_id" % i] = models.PositiveIntegerField()
    seq_dict["r%i_index" % i] = models.PositiveIntegerField()
    seq_dict["r%i_newID" % i] = models.IntegerField()
    seq_dict["r%i_oldID" % i] = models.CharField(max_length=5)
    seq_dict["r%i_ss" % i] = models.CharField(max_length=1, choices=SS_CHOICES)
    seq_dict["r%i_terminal_flag" % i] = models.BooleanField()
    seq_dict["r%i_xpr" % i] = models.BooleanField() # probably should be replaced

    # the loops here are just to save on space/typing
    for j in range(1,8):
        seq_dict["r%i_a%i" % (i,j)] = models.FloatField()
    for j in range(1,6):
        seq_dict["r%i_L%i" % (i,j)] = models.FloatField()
    for j in ("phi", "psi", "ome", "chi", "bm", "bs", "bg", "h_bond_energy", "zeta"):
        seq_dict["r%i_%s" % (i,j)] = models.FloatField()

# Create the Segment model with the fields from the dict
Segment = type('Segment', (Segment_abstract,), seq_dict)
