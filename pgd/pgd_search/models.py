from django.db import models
from django.contrib.auth.models import User
from pgd_core.models import Protein,Residue
from constants import AA_CHOICES, SS_CHOICES, Subscripter
from exceptions import AttributeError
import re
import math

import dbsettings
from dbsettings.loading import set_setting_value

# Search settngs
class SearchSettings(dbsettings.Group):
    segmentSize          = dbsettings.IntegerValue('Current Segment Size', 'Maximum size for segment searches')
    requestedSegmentSize = dbsettings.IntegerValue('Requested Segment Size','Requested size for segment searches.  This value is used to generate tables and data prior to a search of this size is available')
searchSettings = SearchSettings('Splicer')

#set defaults for settings
if not searchSettings.segmentSize:
    set_setting_value('pgd_search.models', '', 'segmentSize', 10)
if not searchSettings.requestedSegmentSize:
    set_setting_value('pgd_search.models', '', 'requestedSegmentSize', 10)


# Search
# A search query submitted by a user
class Search(models.Model):
    user             = models.ForeignKey(User, null=True)
    codes_include    = models.BooleanField(null=True)
    _querySet = None

    # returns the query set that represents this search
    def querySet(self):
        #create querySet if not needed
        if not self._querySet:
            #_querySet = parse_search(self)

            #TODO until the _queryParser is working and merged in just return the first 50 full length segments
            self._querySet = Segment.objects.all().filter(length=10)[:500]

        return self._querySet

# Search_code
# Codes for the proteins searched on
class Search_code(models.Model):
    search  = models.ForeignKey(Search, related_name='codes')
    code    = models.CharField(max_length=4)

# Search_residue
# The search fields per residue
class Search_residue(models.Model):
    search          = models.ForeignKey(Search, related_name='residues')
    index           = models.IntegerField()
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

    # '<field>_include' boolean determines how its query field should be handled
    # Null  : field not included
    # True  : field included as a positive assertion
    # False : field included as a negative assertion
    aa_int_include          = models.BooleanField(null=True)
    a1_include              = models.BooleanField(null=True)
    a2_include              = models.BooleanField(null=True)
    a3_include              = models.BooleanField(null=True)
    a4_include              = models.BooleanField(null=True)
    a5_include              = models.BooleanField(null=True)
    a6_include              = models.BooleanField(null=True)
    a7_include              = models.BooleanField(null=True)
    L1_include              = models.BooleanField(null=True)
    L2_include              = models.BooleanField(null=True)
    L3_include              = models.BooleanField(null=True)
    L4_include              = models.BooleanField(null=True)
    L5_include              = models.BooleanField(null=True)
    ss_include              = models.BooleanField(null=True)
    phi_include             = models.BooleanField(null=True)
    psi_include             = models.BooleanField(null=True)
    ome_include             = models.BooleanField(null=True)
    chi_include             = models.BooleanField(null=True)
    bm_include              = models.BooleanField(null=True)
    bs_include              = models.BooleanField(null=True)
    bg_include              = models.BooleanField(null=True)
    h_bond_energy_include   = models.BooleanField(null=True)
    zeta_include            = models.BooleanField(null=True)
    terminal_flag_include   = models.BooleanField(null=True)


    def __init__(self):
        models.Model.__init__(self)

        # populate 'aa' with a dictionary of allowed values from AA_CHOICES
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
