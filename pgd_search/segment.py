from django.db import models

import math
import re

from pgd_constants import SS_CHOICES, AA_CHOICES
from pgd_core.models import *
from models import searchSettings

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

        except TypeError: # it wasn't an int, it must be a slice
            #return the range of segments specified.  retrieve items through __getitem__() 
            #so items are initialized if needed.
            return [self.__getitem__(z) for z in range(i.start, i.stop, i.step if i.step else 1) ]


    def __setitem__(self, i, v):
        #populate dict and set the FK
        if v:
            self.parent.__dict__['r%i' % i] = v
            self.parent.__dict__['r%i_id' % i] = v.id

        #v is None, clear the object from the dict and clear the FK (if not already None)
        elif self.parent.__dict__.has_key('r%i' % i):
            del self.parent.__dict__['r%i' % i]
            self.parent.__dict__['r%i_id' % i] = None

    def __iter__(self):
        # This function makes a generator object
        def residue_generator(outer):
            for i in range(searchSettings.segmentSize):
                try: # Get the next object...
                    yield outer.__getitem__(i)
                except IndexError: # ...until there are no more.
                    raise StopIteration
        return residue_generator(self)


#this pattern matches any property that is proxied to a residue
proxyPattern = re.compile('^r([\d]+)_(?!id$)([\w]+)')


# Segment_abstract
# A base class for the Segment object
iIndex = int(math.ceil(searchSettings.segmentSize/2.0)-1)
class Segment_abstract(models.Model):

    protein = models.ForeignKey(Protein)
    chainID = models.CharField(max_length=1)
    length  = models.PositiveIntegerField()

    def __init__(self, *args, **kwargs):
        super(Segment_abstract, self).__init__(*args, **kwargs)
        Residue_subscripter('residues', self)

    def __getattribute__(self,name):
        if name == 'offset':
            return 0-self.__dict__['r%i_chainIndex' % iIndex]
    
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
for i in range(searchSettings.segmentSize):

    seq_dict["r%i_id" % i]              = models.PositiveIntegerField(null=True)
    seq_dict["r%i_chainIndex" % i]      = models.PositiveIntegerField(null=True)
    seq_dict["r%i_ss" % i]              = models.CharField(max_length=1, choices=SS_CHOICES, null=True)
    seq_dict["r%i_aa" % i]              = models.CharField(max_length=1, choices=AA_CHOICES, null=True)
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