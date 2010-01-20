from django.db import models
from django.contrib.auth.models import User
from pgd_core.models import Protein,Residue
from pgd_constants import AA_CHOICES, SS_CHOICES, Subscripter
from exceptions import AttributeError
from django import forms
import re
from math import ceil

import dbsettings

from django.db.models import Q

range_re = re.compile("(?<=[^-<>=])-")
comp_re  = re.compile("^([<>]=?)?")


class SearchSettings(dbsettings.Group):
    """
    Search Settings
    """
    query_limit          = dbsettings.IntegerValue('Result Set Limit', 'Maximum size of results.  Due to the heavy cpu requirement for calculating statistics this option limits how many results can be returned by a search', default=50000000)
    segmentSize          = dbsettings.IntegerValue('Current Segment Size', 'Maximum size for segment searches', default=10)
    requestedSegmentSize = dbsettings.IntegerValue('Requested Segment Size','Requested size for segment searches.  This value is used to generate tables and data prior to a search of this size is available', default=10)
searchSettings = SearchSettings('Search')


class Search(models.Model):
    """
    A query submitted by a user
    """
    dataset_version  = models.CharField(max_length='100')
    isPublic         = models.BooleanField()
    timestamp        = models.DateTimeField(null=True)
    title            = models.CharField(max_length='300')
    description      = models.CharField(max_length='5000')
    user             = models.ForeignKey(User, null=True)
    codes_include    = models.NullBooleanField(null=True)
    threshold        = models.IntegerField(null=True)
    resolution_min   = models.FloatField(null=True)
    resolution_max   = models.FloatField(null=True)
    rfactor_min      = models.FloatField(null=True)
    rfactor_max      = models.FloatField(null=True)
    rfree_min        = models.FloatField(null=True)
    rfree_max        = models.FloatField(null=True)

    segmentLength    = models.IntegerField()
    _querySet = None

    def querySet(self):
        """
        returns the query set that represents this search 
        """

        #create querySet if not needed
        #if not self._querySet:
        #    self._querySet = self.parse_search()

        #return self._querySet

        # for now parse query every time.  otherwise the query results will be stored in the session
        return self.parse_search()

    # Return the segments matched by the search.
    def parse_search(self):

        # Start with all segments...
        query = Residue.objects.all()

        # ...filter by segmentLength...
        #if self.segmentLength > 1:
        #    query = query.filter(length__gte=self.segmentLength)

        # ...filter by code lists...
        if self.codes_include != None:
            if self.codes_include:
                query = query.filter(protein__in=(x.code for x in self.codes.all()))
            else:
                for x in self.codes.all():
                    query = query.exclude(protein=x.code)

        # ...filter by code lists...
        if self.resolution_min != None:
            query = query.filter(protein__resolution__gte=self.resolution_min)

        # ...filter by resolution...
        if self.resolution_max != None:
            query = query.filter(protein__resolution__lte=self.resolution_max)

        # ...filter by rfactor...
        if self.rfactor_min != None:
            query = query.filter(protein__rfactor__gte=self.rfactor_min)

        # ...filter by rfactor...
        if self.rfactor_max != None:
            query = query.filter(protein__rfactor__lte=self.rfactor_max)

        #...filter by rfree...
        if self.rfree_min != None:
            query = query.filter(protein__rfree__gte=self.rfree_min)

        # ...filter by rfree...
        if self.rfree_max != None:
            query = query.filter(protein__rfree__lte=self.rfree_max)


        # ...filter by threshold...
        if self.threshold != None:
            query = query.filter(protein__threshold__lte=self.threshold)

        # ...filter by query strings (for values and value ranges)...
        def compare(x,y):
            if x == y:
                return 0
            elif x < y:
                return -1
            return 1
        residues = sorted(self.residues.all(), compare, lambda x: abs(x.index) )

        for search_res in residues: # iterate through all search residues in self
            #seg_prefix = "r%i_"%(
            #    # convert from the search residue index to indexes 0...n
            #    search_res.index+int(ceil(searchSettings.segmentSize/2.0)-1)
            #)

            # get field prefix for this residue
            if search_res.index == 0:
                seg_prefix = ''
            elif search_res.index < 0:
                seg_prefix = ''.join(['prev__' for i in range(search_res.index, 0)])
            else:
                seg_prefix = ''.join(['next__' for i in range(search_res.index)])

            # ...handle boolean values...
            #   (Filter on each binary field that is not set to NULL/None.)
            query = query.filter(**dict((
                    (seg_prefix+field, search_res.__dict__[field])
                    for field in (
                        'xpr',
                    ) if search_res.__dict__[field] != None
                )))

            # ...handle the '_int' values...
            #   ('_int' values are a series of booleans stored grouped in an integer.)
            for field,choices in filter(
                   # use only those (_int,_CHOICES) pairs with a '_include' value.
                   lambda x: search_res.__dict__[x[0]+'_include'] != None,
                   (
                       ('aa_int', AA_CHOICES),
                       ('ss_int', SS_CHOICES),
                   )
                ):
                query = query.__getattribute__(
                    # call either 'filter' or 'exclude', depending on the value of '_include'
                    'filter' if search_res.__dict__[field+'_include'] else 'exclude'
                )(
                    # check to see that the value of the segment residue is in the set of
                    # residues described in the '_int' of the search residue.
                    **{seg_prefix+field[0:-4]+"__in": [
                        # get the designated choice names as stored in the database.
                        choice[0] for index,choice in enumerate(choices) if search_res.__dict__[field]&1<<index
                    ]}
                )

            # ...handle query strings...
            for field in filter(
                # use only the fields with a '_include' value. 
                lambda x: search_res.__dict__[x+'_include'] != None,
                (
                    'a1',   'a2',   'a3',   'a4',   'a5',   'a6',   'a7',
                    'L1',   'L2',   'L3',   'L4',   'L5',
                    'phi',  'psi',  'ome',  'chi1', 'chi2', 'chi3', 'chi4',
                    'bm',   'bs',   'bg',
                    'h_bond_energy',
                    'zeta',
                )
            ):

                # seg_field is the name of the property of the given residue in the database
                seg_field = seg_prefix+field

                constraints = []

                for constraint in str(search_res.__dict__[field]).split(','):

                    # The constraint may be a range...
                    if range_re.search(constraint):

                        min,max = [float(lim) for lim in range_re.split(constraint)]

                        limits = (
                            Q(**{seg_field+'__gte' : float(min)}),
                            Q(**{seg_field+'__lte' : float(max)}),
                        )

                        constraints.append(
                            # Apply 'or' logic for a wraparound range
                            # or apply 'and' logic for a regular range
                            (limits[0] & limits[1]) if (min <= max) else (limits[0] | limits[1])
                        )
                    # ...or the constraint may be a value comparison.
                    else:

                        constraints.append(Q(**{
                            seg_field + {
                                ''   : '',
                                '>'  : '__gt',
                                '>=' : '__gte',
                                '<'  : '__lt',
                                '<=' : '__lte',
                            }[comp_re.search(constraint).group(0)] :
                            # extract the numeric value
                            constraint.strip('><=')
                        }))

                query = query.__getattribute__(
                    # call either 'filter' or 'exclude', depending on the value of '_include'
                    'filter' if search_res.__dict__[field+'_include'] else 'exclude'
                )(
                    # OR each of the statements in the query string together
                    reduce(
                        lambda x,y: x|y,
                        constraints
                    )
                )
        return query


class Search_code(models.Model):
    """
    Codes indicating to which proteins a Search is applied.
    """
    search  = models.ForeignKey(Search, related_name='codes')
    code    = models.CharField(max_length=4)


class Search_residue(models.Model):
    """
    The per-residue properties of a Search
    """
    search          = models.ForeignKey(Search, related_name='residues')
    index           = models.IntegerField()
    chainID         = models.CharField(max_length=1, null=True)

    # Performing bitwise operations on aa_int gives the set of amino acids to check against.
    # aa_int&1<<i   AA_CHOICE[i]
    # 0 : Not in set
    # 1 : In set
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

    # Performing bitwise operations on ss_int gives the set of secondary structures to check against.
    # ss_int&1<<i   AA_CHOICE[i]
    # 0 : Not in set
    # 1 : In set
    ss_int          = models.IntegerField(null=True)

    phi             = models.CharField(max_length=30, null=True)
    psi             = models.CharField(max_length=30, null=True)
    ome             = models.CharField(max_length=30, null=True)
    chi1            = models.CharField(max_length=30, null=True)
    chi2            = models.CharField(max_length=30, null=True)
    chi3            = models.CharField(max_length=30, null=True)
    chi4            = models.CharField(max_length=30, null=True)
    bm              = models.CharField(max_length=30, null=True)
    bs              = models.CharField(max_length=30, null=True)
    bg              = models.CharField(max_length=30, null=True)
    h_bond_energy   = models.CharField(max_length=30, null=True)
    zeta            = models.CharField(max_length=30, null=True)
    xpr             = models.NullBooleanField(null=True) # this field may not be necessary; it has never been implemented

    # '<field>_include' boolean determines how its query field should be handled
    # Null  : field not included
    # True  : field included as a positive assertion
    # False : field included as a negative assertion
    aa_int_include          = models.NullBooleanField(null=True)
    a1_include              = models.NullBooleanField(null=True)
    a2_include              = models.NullBooleanField(null=True)
    a3_include              = models.NullBooleanField(null=True)
    a4_include              = models.NullBooleanField(null=True)
    a5_include              = models.NullBooleanField(null=True)
    a6_include              = models.NullBooleanField(null=True)
    a7_include              = models.NullBooleanField(null=True)
    L1_include              = models.NullBooleanField(null=True)
    L2_include              = models.NullBooleanField(null=True)
    L3_include              = models.NullBooleanField(null=True)
    L4_include              = models.NullBooleanField(null=True)
    L5_include              = models.NullBooleanField(null=True)
    ss_int_include          = models.NullBooleanField(null=True)
    phi_include             = models.NullBooleanField(null=True)
    psi_include             = models.NullBooleanField(null=True)
    ome_include             = models.NullBooleanField(null=True)
    chi1_include             = models.NullBooleanField(null=True)
    chi2_include             = models.NullBooleanField(null=True)
    chi3_include             = models.NullBooleanField(null=True)
    chi4_include             = models.NullBooleanField(null=True)
    bm_include              = models.NullBooleanField(null=True)
    bs_include              = models.NullBooleanField(null=True)
    bg_include              = models.NullBooleanField(null=True)
    h_bond_energy_include   = models.NullBooleanField(null=True)
    zeta_include            = models.NullBooleanField(null=True)


    def __init__(self, *args, **kwargs):
        models.Model.__init__(self, *args, **kwargs)

        # populate 'aa' with a dictionary of allowed values from AA_CHOICES
        self.aa = dict([(aa_choice[1],1 if self.aa_int == None else 1&self.aa_int>>aa_index) for aa_index,aa_choice in enumerate(AA_CHOICES)])
        # populate 'ss' with a dictionary of allowed values from SS_CHOICES
        self.ss = dict([(ss_choice[1],1 if self.ss_int == None else 1&self.ss_int>>ss_index) for ss_index,ss_choice in enumerate(SS_CHOICES)])



class ResidueProxy():
    """
    This is a proxy to properties stored in the segment.  It is used
    to emulate having a Residue but really just returns properties
    that are stored in the segment object.  This is used because
    it is faster than fetching the Residue from the database

    This class is used in conjunction with ResidueProxy_subscripter
    """
    def __init__(self, index, parent):
        self.parent = parent
        self.string = 'r%i_%%s' % index

    def __getitem__(self, i):
        return self.parent.__dict__[self.string % i]


class ResidueProxy_subscripter():
    """
    ResidueProxy_subscripter

    This class emulates a list of residues.  It returns ResidueProxy
    objects that emulate a real residue object.
    """

    def __init__(self, key, parent):
        self.parent = parent
        #add this instance to the parent. doing this here
        #makes defining subscriptor instance simpler because
        #you only need to specify the key once
        parent.__dict__[key] = self

        # create a list of proxy objects for the residues in this segment
        self.proxies = []
        l_append = self.proxies.append
        for i in range(searchSettings.segmentSize):
            l_append(ResidueProxy(i, parent))

    def __getitem__(self, i):
        try: # Get the object, if it's been instantiated...
            return self.proxies[i]

        except TypeError: # it wasn't an int, it must be a slice
            #return the range of segments specified.  retrieve items through __getitem__() 
            #so items are initialized if needed.
            return self.proxies[i.start:i.stop:i.step if i.step else 1]


class Residue_subscripter():
    """
    This is a lazy loading list of residues that a segment contains
    When residues are requested the id's are fetched from the parent
    segment and used to lookup the residue.  Residues are cached so
    repeat lookups do not require database interaction

    additionally this class also allows iteration of residues
    """

    def __init__(self, key, parent):
        self.parent = parent
        #add this instance to the parent. doing this here
        #makes defining subscriptor instance simpler because
        #you only need to specify the key once
        parent.__dict__[key] = self


    def __getitem__(self, i):
        """
        called when an index is requested (ie. residues[1])
        """
        try: # Get the object, if it's been instantiated...
            return self.parent.__dict__['r%i' % i]
            #return self.foo[i]

        except KeyError: # ...otherwise, instantiate it.
            pass
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
            return self.foo[i.start:i.stop:i.step if i.step else 1]


    def __setitem__(self, i, v):
        """
        called when an index is set (ie. residues[1]=foo)
        """

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


class SegmentedProtein(models.Model):
    """
    Parallel to the Protein model to add fields specific to the segmentation
    of the protein chain into segment objects

    This class is not an extension of the Protein model because it is intended
    to have limited use.  Linking this to the Protein object creates issues
    creating and querying (extra joins) Protein and SegmentedProtein objects
    """
    code = models.CharField(max_length=4, primary_key=True, unique=True)

    # Date of the pdb the segments were generated from.  This is different from
    # Protein.pdb_date.  It serves the same purpose of indicating if a Segment
    # is up to date.  The difference is that this tracks updates for the
    # Segment table.  The import process has separate transactions for the core
    # protein and segments.  Its possible for the two to get out of sync if an
    # import fails while proteins are being processed.
    pdb_date = models.DateTimeField()


""" ====================================================================
Segment Model

This model consists of 2 parts:
   1) an abstract base class that contains all of the normal properties
       and function definitions
   2) a dynamically created child class that contains a set of fields
      that dynamically resize depending on maximum possible segment
      length.  This works by defining a dictionary of fields and 
      passing it to the type() function which creates a class bound
      to the namespace of this file

This allows the segment model to change size dynamically at runtime.

==================================================================== """
#this pattern matches any property that is proxied to a residue
proxyPattern = re.compile('^r([\d]+)_(?!id$)([\w]+)')


iIndex = int(ceil(searchSettings.segmentSize/2.0)-1)
class Segment_abstract(models.Model):

    protein = models.ForeignKey(Protein)
    chainID = models.CharField(max_length=1)
    length  = models.PositiveIntegerField()

    def __init__(self, *args, **kwargs):
        super(Segment_abstract, self).__init__(*args, **kwargs)
        Residue_subscripter('residues', self)
        ResidueProxy_subscripter('residueProxies', self)

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
    seq_dict["r%i_oldID" % i]           = models.CharField(max_length=5, null=True)
    seq_dict["r%i_chainIndex" % i]      = models.PositiveIntegerField(null=True)
    seq_dict["r%i_ss" % i]              = models.CharField(max_length=1, choices=SS_CHOICES, null=True)
    seq_dict["r%i_aa" % i]              = models.CharField(max_length=1, choices=AA_CHOICES, null=True)
    seq_dict["r%i_xpr" % i]             = models.NullBooleanField(null=True) # probably should be replaced

    # the loops here are just to save on space/typing
    for j in range(1,8):
        seq_dict["r%i_a%i" % (i,j)] = models.FloatField(null=True)
    for j in range(1,6):
        seq_dict["r%i_L%i" % (i,j)] = models.FloatField(null=True)
    for j in ("phi", "psi", "ome", "chi", "bm", "bs", "bg", "h_bond_energy", "zeta"):
        seq_dict["r%i_%s" % (i,j)] = models.FloatField(null=True)

# Create the Segment model class with the fields from the dict
Segment = type('Segment', (Segment_abstract,), seq_dict)

class saveSearchForm(forms.Form):
    title       = forms.CharField(label='Title')
    description = forms.CharField(label='Description', widget=forms.Textarea)
    isPublic    = forms.BooleanField(label='Publically Viewable',required=False)
    id          = forms.IntegerField(None, widget=forms.HiddenInput, required=False)

