from exceptions import AttributeError
from math import ceil
import re
import cPickle

from django import forms
from django.conf import settings
from django.db import models
from django.db.models import Q
from django.contrib.auth.models import User

from pgd_core.models import Protein,Residue
from pgd_constants import AA_CHOICES, AA_CHOICES_DICT, SS_CHOICES
from pgd_splicer.sidechain import bond_lengths_string_dict, bond_angles_string_dict
from pgd_core.util import residue_indexes


range_re = re.compile("(?<=[^-<>=])-")
comp_re  = re.compile("^([<>]=?)?")


class RDict(dict):
    """ Helper class for accessing a dict's properties as if they were member
    variables.  also return None for any value not found as a member
    """
    def __init__(self, dict_):
        self.update(dict_)
        super(dict, self).__init__()

    def __getitem__(self, k):
        try:
            return super(RDict, self).__getitem__(k)
        except KeyError:
            return None

    def __getattribute__(self, k):
        try:
            return super(dict, self).__getattribute__(k)
        except AttributeError:
            try:
                return self[k]
            except KeyError:
                return None


class Search(models.Model):
    """
    A query submitted by a user
    """
    dataset_version  = models.CharField(max_length='100')
    isPublic         = models.BooleanField(default=True)
    timestamp        = models.DateTimeField(null=True)
    title            = models.CharField(max_length='300')
    description      = models.CharField(max_length='5000')
    user             = models.ForeignKey(User, null=True)
    data_internal    = models.TextField()
    __data = None

    @property
    def data(self):
        if self.__data is None and self.data_internal:
            self.__data = cPickle.loads(str(self.data_internal))
        return self.__data

    @data.setter
    def data(self, value):
        self.__data = value
        self.data_internal = None

    @property
    def segmentLength(self):
        return int(self.data['residues'])


    @property
    def residues(self):
        """ iterates residues found in the search.  The objects yielded are
        self.data wrapped in a Segmenter.
        """
        data = self.data
        for i in range(self.segmentLength):
            yield Segmenter(data, i)



    def save(self):
        """ serialize parameters pre-save, object probably wont't be saved often """
        if self.data_internal is None and self.__data is not None:
            self.data_internal = cPickle.dumps(self.__data)
        super(Search, self).save()

    _querySet = None

    def querySet(self):
        """
        returns the query set that represents this search
        """
        # for now parse query every time.  otherwise the query results will be stored in the session
        return self.parse_search()

    # Return the segments matched by the search.
    def parse_search(self):

        # Start with all segments...
        query = Residue.objects.all()
        data = self.data
        if not data:
            # if no params return everything
            return query
        data = RDict(data)

        # ...filter by code lists...
        if data.proteins_i != None:
            #codes = data.proteins.replace(' ','').rstrip(',').split(',')
            data.proteins = data.proteins.replace(' ','').replace(',', '')
            codes = [data.proteins[i:i+4] for i in range(0, len(data.proteins), 4)]
            if data.proteins_i:
                query = query.filter(protein__code__in=codes)
            else:
                query = query.exclude(protein__code__in=codes)

        # ...filter by code lists...
        if data.resolutionMin != None:
            query = query.filter(protein__resolution__gte=data.resolutionMin)

        # ...filter by resolution...
        if data.resolutionMax != None:
            query = query.filter(protein__resolution__lte=data.resolutionMax)

        # ...filter by rfactor...
        if data.rfactorMin != None:
            query = query.filter(protein__rfactor__gte=data.rfactorMin)

        # ...filter by rfactor...
        if data.rfactorMax != None:
            query = query.filter(protein__rfactor__lte=data.rfactorMax)

        #...filter by rfree...
        if data.rfreeMin != None:
            query = query.filter(protein__rfree__gte=data.rfreeMin)

        # ...filter by rfree...
        if data.rfreeMax != None:
            query = query.filter(protein__rfree__lte=data.rfreeMax)


        # ...filter by threshold...
        if data.threshold != None:
            query = query.filter(protein__threshold__lte=data.threshold)

        # ...filter by query strings (for values and value ranges)...
        def compare(x,y):
            if x == y:
                return 0
            elif x < y:
                return -1
            return 1
        #residues = sorted(data.residues, compare, lambda x: abs(x.index) )


        indexes = residue_indexes(int(data.residues))
        for index in indexes: # iterate through all search residues in self
            search_res = Segmenter(data, index)

            # get field prefix for this residue
            if index == 0:
                seg_prefix = ''
            elif index < 0:
                seg_prefix = ''.join(['prev__' for i in range(index, 0)])
            else:
                seg_prefix = ''.join(['next__' for i in range(index)])

            # add isnull to ensure segment length is correct.  segments missing
            # a residue would be too short
            if seg_prefix != '':
                query = query.filter(**{'%sisnull' % seg_prefix:False})

            # get possible sidechain fields based on selected AA types
            sidechain_fields = []
            """for aa_type in filter(lambda x: x[1],search_res.aa.items()):
                aa_type_upper = aa_type[0].upper()
                for field in bond_lengths_string_dict[aa_type_upper]:
                    sidechain_fields.append('sidechain_%s__%s' % (aa_type_upper, field))
                for field in bond_angles_string_dict[aa_type_upper]:
                    sidechain_fields.append('sidechain_%s__%s' % (aa_type_upper, field))
            """

            # ...handle boolean values...
            #   (Filter on each binary field that is not set to NULL/None.)
            # XXX there don't appear to be any of these right now, xpr isn't on the form
            '''query = query.filter(**dict((
                    (seg_prefix+field, search_res.__dict__[field])
                    for field in (
                        'xpr',
                    ) if search_res.__dict__[field] != None
                )))
            '''

            # ...handle the '_int' values...
            #   ('_int' values are a series of booleans stored grouped in an integer.)
            for field,choices in filter(
                   # use only those (_int,_CHOICES) pairs with a '_include' value.
                   lambda x: search_res.__getitem__(x[0]+'_i') != None,
                   (
                       ('aa', AA_CHOICES),
                       ('ss', SS_CHOICES),
                   )
                ):
                query = query.__getattribute__(
                    # call either 'filter' or 'exclude', depending on the value of '_include'
                    'filter' if search_res.__getitem__(field+'_i') else 'exclude'
                )(
                    # check to see that the value of the segment residue is in the set of
                    # residues described in the '_int' of the search residue.
                    **{"%s%s__in" % (seg_prefix,field): search_res.__getitem__(field)}
                )

            # ...handle query strings...
            fields = filter(
                # use only the fields with a '_include' value.
                lambda x: search_res.__getitem__(x+'_i') != None,
                (
                    'a1',   'a2',   'a3',   'a4',   'a5',   'a6',   'a7',
                    'L1',   'L2',   'L3',   'L4',   'L5',
                    'phi',  'psi',  'ome',  'chi1', 'chi2', 'chi3', 'chi4', 'chi5',
                    'bm',   'bs',   'bg',   'occm', 'occscs',
                    'h_bond_energy',
                    'zeta',
                )
            )
            query = self.filter_fields(fields, query, search_res, seg_prefix)

            # ... handle sidechain query strings ...
            sidechain_fields = []
            if search_res.aa:
                for aa_type in [AA_CHOICES_DICT[aa].upper() for aa in search_res.aa]:
                    if aa_type in bond_lengths_string_dict:
                        field_base = '%s__%%s' % aa_type
                        for field in bond_lengths_string_dict[aa_type]:
                            key = field_base % field
                            if search_res[key]:
                                sidechain_fields.append(key)

                        for field in bond_angles_string_dict[aa_type]:
                            key = field_base % field
                            if search_res[key]:
                                sidechain_fields.append(key)

            seg_prefix = '%ssidechain_' % seg_prefix
            query = self.filter_fields(sidechain_fields, query, search_res, seg_prefix)

        return query

    def filter_fields(self, fields, query, search_res, seg_prefix):
        """
        Filters the fields passed in
        """

        for field in fields:
            # seg_field is the name of the property of the given residue in the database
            seg_field = seg_prefix+field

            constraints = []

            for constraint in str(search_res.__getitem__(field)).split(','):

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
                'filter' if search_res.__getitem__(field+'_i') else 'exclude'
            )(
                # OR each of the statements in the query string together
                reduce(
                    lambda x,y: x|y,
                    constraints
                )
            )

        return query


class Segmenter(object):
    """ segments form data for easier access """
    def __init__(self,dict__, i):
        self.i = i
        self.__dict = dict__

    def __contains__(self, k):
        return '%s_%d' % (k, self.i) in self.__dict

    def __getitem__(self, k):
        return self.__dict['%s_%d' % (k, self.i)]

    def __getattribute__(self, k):
        try:
            return super(Segmenter, self).__getattribute__(k)
        except AttributeError:
            return self.__dict['%s_%d' % (k, self.i)]

    def __str__(self):
        return str(self.__dict)


class Search_code(models.Model):
    """
    Codes indicating to which proteins a Search is applied.
    """
    search  = models.ForeignKey(Search, related_name='codes')
    code    = models.CharField(max_length=4)


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
        for i in range(settings.SEGMENT_SIZE):
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
            for i in range(settings.SEGMENT_SIZE):
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


iIndex = int(ceil(settings.SEGMENT_SIZE/2.0)-1)
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
for i in range(settings.SEGMENT_SIZE):

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
