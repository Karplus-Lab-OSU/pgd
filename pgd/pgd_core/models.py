import math

from django.db import models
from pgd_constants import AA_CHOICES, SS_CHOICES, AA_CHOICES_DICT

# Protein model
# (was 'protein_info')
# Contains the general information about a protein
class Protein(models.Model):
    code        = models.CharField(max_length=4, primary_key=True, unique=True)
    threshold   = models.IntegerField() # new type; this should probably be a boolean type
    resolution  = models.FloatField(db_index=True)
    rfactor     = models.FloatField(db_index=True)
    rfree       = models.FloatField(db_index=True)

    # date the source PDB file was created.  This property is used to check for
    # updates allowing up-to-date proteins to be skipped.
    pdb_date    = models.DateTimeField()

    def __unicode__(self):
        return self.code

# Chain model
# Contains information about chains within a protein
# 
# NOTE: This model is in a temporary state that allows easier conversion from the current database
# the primary key will be a 5 character string derived from the chainID and protein code.
# this will all a very simple query to generate both the table of chains and set the FK
# in the residue table.  Once splicer is updated and importing directly to the database the PK
# will be changed to an integer.
class Chain (models.Model):
        id      = models.CharField(max_length=5, primary_key=True)
        protein = models.ForeignKey(Protein, related_name='chains')
        code    = models.CharField(max_length=1)


# (Note: fields need to be commented)
class Residue(models.Model):
    """
    Residue model - Contains information regarding each amino acid and its geometry

    Indexing:
     -chainIndex:    An integer based identifier scheme for residues in the
                     chain.  The residues are numbered started with 1 and skip
                     a number for any chain breaks.  This is used internally
                     for querying residues that are next to eachother in the
                     chain.

     -oldID:          Identifier taken from the PDB file.  This is a composite
                     field taken from the residue_id and insertion code if any.
                     This field is displayed to the user because it corresponds
                     to the identifier in the PDB which they will need to do
                     further research on a protein.

     -terminal_flag: A flag indicating a residue is next to a chain break.
                     This flag makes it possible to quickly search for or
                     identify a chain break without comparing the next residue
    """

    protein         = models.ForeignKey(Protein, related_name='residues')
    chain           = models.ForeignKey(Chain, related_name='residues')
    prev            = models.ForeignKey('self', related_name='prev_next')
    next            = models.ForeignKey('self', related_name='next_prev')
    aa              = models.CharField(max_length=1, choices=AA_CHOICES) # new type
    chainID         = models.CharField(max_length=1) # integer id
    oldID           = models.CharField(max_length=5, null=True)# id[icode] from pdb file
    chainIndex      = models.PositiveIntegerField()
    a1              = models.FloatField(null=True)
    a2              = models.FloatField()
    a3              = models.FloatField()
    a4              = models.FloatField()
    a5              = models.FloatField()
    a6              = models.FloatField(null=True)
    a7              = models.FloatField(null=True)
    L1              = models.FloatField(null=True)
    L2              = models.FloatField()
    L3              = models.FloatField()
    L4              = models.FloatField()
    L5              = models.FloatField()
    ss              = models.CharField(max_length=1, choices=SS_CHOICES) # new type (was blob, but all entries 1 char)
    phi             = models.FloatField()
    psi             = models.FloatField()
    ome             = models.FloatField()
    chi             = models.FloatField()
    bm              = models.FloatField()
    bs              = models.FloatField()
    bg              = models.FloatField(null=True)
    h_bond_energy   = models.FloatField()
    zeta            = models.FloatField()
    terminal_flag   = models.BooleanField(default=False)#indicates this residue is next to a chain break
    xpr             = models.BooleanField() # this field may not be necessary; it has never been implemented

    def __init__(self, *args, **kwargs):
        self.segment = Segmenter(self)
        models.Model.__init__(self, *args, **kwargs)

    #def __str__(self):
    #    return '%d' % self.chainIndex

    def __getattribute__(self,name):
        if name == 'aa_full':
            return AA_CHOICES_DICT[self.aa]
        # normal attribute
        else:
            return object.__getattribute__(self, name)


class Segmenter():
    """
    Helper Class for walking a protein chain through prev/next properties
    of Residue
    """
    def __init__(self, residue):
        self.residue = residue

    def __getitem__(self, index):
        if not type(index) == int:
            raise IndexError
        if index == 0:
            return self.residue
        residue = self.residue
        try:
            if index < 0:
                while index != 0 and residue:
                    residue = residue.prev
                    index += 1
            else:
                while index != 0 and residue:
                    residue = residue.next
                    index -= 1
        except Residue.DoesNotExist:
            return None
        return residue


def determine_alias(query, index):
        """
        determines the table alias used for a given residue index.
        
        XXX This takes into account django internal structure as of 12/29/2009
        this may change with future releases.
        
        query.join_map is a dict mapping a tuple of (table1, table2, fk, key)
        mapped to a list of aliases the table is joined on.  multiple aliases
        means the table was joined on itself multiple times.
        
        we must walk the list of joins to find the index number we want.
        
        @returns alias if table is joined, otherwise None
        """
        query = query.query
        if index == 0:
            return 'pgd_core_residue'
        if index > 0:
            k = ('pgd_core_residue','pgd_core_residue','next_id','id')
        else:
            k = ('pgd_core_residue','pgd_core_residue','prev_id','id')
            
        if not query.join_map.has_key(k):
            return None
        try:
            return query.join_map[k][int(math.fabs(index))-1]
        except IndexError:
            return None
