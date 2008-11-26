from django.db import models
from django.contrib.auth.models import User
from pgd_core.models import Protein,Residue
from constants import AA_CHOICES, SS_CHOICES, SEQUENCE_SIZE, Subscripter
from exceptions import AttributeError
import re

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
    chainID         = models.CharField(max_length=1, null=True)
    aa_int          = models.IntegerField(null=True)

    a1              = models.CharField(max_length=30, null=True)
    a2              = models.CharField(max_length=30, null=True)
    a3              = models.CharField(max_length=30, null=True)
    a4              = models.CharField(max_length=30, null=True)
    a5              = models.CharField(max_length=30, null=True)
    a6              = models.CharField(max_length=30, null=True)
    a7              = models.CharField(max_length=30, null=True)
    L1              = models.CharField(max_length=30, null=True)
    L2              = models.CharField(max_length=30, null=True)
    L3              = models.CharField(max_length=30, null=True)
    L4              = models.CharField(max_length=30, null=True)
    L5              = models.CharField(max_length=30, null=True)
    ss              = models.CharField(max_length=1, choices=SS_CHOICES, null=True) # new type (was blob, but all entries 1 char)
    phi             = models.CharField(max_length=30, null=True)
    psi             = models.CharField(max_length=30, null=True)
    ome             = models.CharField(max_length=30, null=True)
    chi             = models.CharField(max_length=30, null=True)
    bm              = models.CharField(max_length=30, null=True)
    bs              = models.CharField(max_length=30, null=True)
    bg              = models.CharField(max_length=30, null=True)
    h_bond_energy   = models.CharField(max_length=30, null=True)
    zeta            = models.CharField(max_length=30, null=True)
    terminal_flag   = models.BooleanField(null=True)
    xpr             = models.BooleanField(null=True) # this field may not be necessary; it has never been implemented

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
            # If the object can be instantiated, return it
            if self.parent.__dict__['r%i_id' % i]:
                self.parent.__dict__['r%i' % i] = Residue.objects.filter(id=self.parent.__dict__['r%i_id' % i])[0]
                return self.parent.__dict__['r%i' % i]

            # otherwise return None.  These objects are proxies to underlying properties 
            # so they must always return a values or None
            else:
                return None

    def __setitem__(self, i, v):
        #populate dict and set the FK
        if v:
            self.parent.__dict__['r%i' % i] = v
            self.parent.__dict__['r%i_id' % i] = v.id

        #v is None, clear the object from the dict and clear the FK
        else:
            del self.parent.__dict__['r%i' % i]
            self.parent.__dict__['r%i_id' % i] = None

    def __iter__(self):
        # This function makes a generator object
        def residue_generator(outer):
            for i in range(SEQUENCE_SIZE):
                try: # Get the next object...
                    yield outer.__getitem__(i)
                except IndexError: # ...until there are no more.
                    raise StopIteration
        return residue_generator(self)


#this pattern matches any property that is proxied to a residue
proxyPattern = re.compile('^r([\d]+)_(?!id$)([\w]+)')


# Segment_abstract
# A base class for the Segment object
class Segment_abstract(models.Model):

    protein = models.ForeignKey(Protein)
    chainID = models.CharField(max_length=1)

    def __init__(self):
        models.Model.__init__(self)
        Residue_subscripter('residues', self)

    def __getattribute__(self,name):
        #check for properties proxied to a residue object
        match = proxyPattern.match(name)
        if match:
            index = int(match.group(1))
            attr  = match.group(2)
            residue = object.__getattribute__(self, 'residues')[index]
            if residue:
                return residue.__dict__[attr]

            #return none if residue is none (transitive)
            return None

        # not a proxied attribute
        else:
            return object.__getattribute__(self, name)

    class Meta:
        abstract = True

# Build a dict for the fields of variable number
seq_dict = {'__module__' : 'pgd_search.models'}
for i in range(10):

    seq_dict["r%i_id" % i]              = models.PositiveIntegerField(null=True)
    seq_dict["r%i_index" % i]           = models.PositiveIntegerField(null=True)
    seq_dict["r%i_ss" % i]              = models.CharField(max_length=1, choices=SS_CHOICES, null=True)
    seq_dict["r%i_terminal_flag" % i]   = models.BooleanField(null=True)
    seq_dict["r%i_xpr" % i]             = models.BooleanField(null=True) # probably should be replaced

    # the loops here are just to save on space/typing
    for j in range(1,8):
        seq_dict["r%i_a%i" % (i,j)] = models.FloatField(null=True)
    for j in range(1,6):
        seq_dict["r%i_L%i" % (i,j)] = models.FloatField(null=True)
    for j in ("phi", "psi", "ome", "chi", "bm", "bs", "bg", "h_bond_energy", "zeta"):
        seq_dict["r%i_%s" % (i,j)] = models.FloatField(null=True)

# Create the Segment model with the fields from the dict
Segment = type('Segment', (Segment_abstract,), seq_dict)
