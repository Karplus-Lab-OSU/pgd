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


class Sidechain_ARG(models.Model):
    CB_CG = models.FloatField(null=True)
    CG_CD = models.FloatField(null=True)
    CD_NE = models.FloatField(null=True)
    NE_CZ = models.FloatField(null=True)
    CZ_NH1 = models.FloatField(null=True)
    CZ_NH2 = models.FloatField(null=True)
    CA_CB_CG = models.FloatField(null=True)
    CB_CG_CD = models.FloatField(null=True)
    CG_CD_NE = models.FloatField(null=True)
    CD_NE_CZ = models.FloatField(null=True)
    NE_CZ_NH1 = models.FloatField(null=True)
    NE_CZ_NH2 = models.FloatField(null=True)
    NH1_CZ_NH2 = models.FloatField(null=True)


class Sidechain_ASN(models.Model):
    CB_CG = models.FloatField(null=True)
    CG_OD1 = models.FloatField(null=True)
    CG_ND2 = models.FloatField(null=True)
    CA_CB_CG = models.FloatField(null=True)
    CB_CG_OD1 = models.FloatField(null=True)
    CB_CG_ND2 = models.FloatField(null=True)
    OD1_CG_ND2 = models.FloatField(null=True)


class Sidechain_ASP(models.Model):
    CB_CG = models.FloatField(null=True)
    CG_OD1 = models.FloatField(null=True)
    CG_OD2 = models.FloatField(null=True)
    CA_CB_CG = models.FloatField(null=True)
    CB_CG_OD1 = models.FloatField(null=True)
    CB_CG_OD2 = models.FloatField(null=True)
    OD1_CG_OD2 = models.FloatField(null=True)


class Sidechain_CYS(models.Model):
    CB_SG = models.FloatField(null=True)
    CA_CB_SG = models.FloatField(null=True)


class Sidechain_GLN(models.Model):
    CB_CG = models.FloatField(null=True)
    CG_CD = models.FloatField(null=True)
    CD_OE1 = models.FloatField(null=True)
    CD_NE2 = models.FloatField(null=True)
    CA_CB_CG = models.FloatField(null=True)
    CB_CG_CD = models.FloatField(null=True)
    CG_CD_OE1 = models.FloatField(null=True)
    CG_CD_NE2 = models.FloatField(null=True)
    OE1_CD_NE2 = models.FloatField(null=True)


class Sidechain_GLU(models.Model):
    CB_CG = models.FloatField(null=True)
    CG_CD = models.FloatField(null=True)
    CD_OE1 = models.FloatField(null=True)
    CD_OE2 = models.FloatField(null=True)
    CA_CB_CG = models.FloatField(null=True)
    CB_CG_CD = models.FloatField(null=True)
    CG_CD_OE1 = models.FloatField(null=True)
    CG_CD_OE2 = models.FloatField(null=True)
    OE1_CD_OE2 = models.FloatField(null=True)


class Sidechain_HIS(models.Model):
    CB_CG = models.FloatField(null=True)
    CG_ND1 = models.FloatField(null=True)
    ND1_CE1 = models.FloatField(null=True)
    CE1_NE2 = models.FloatField(null=True)
    NE2_CD2 = models.FloatField(null=True)
    CG_CD2 = models.FloatField(null=True)
    CA_CB_CG = models.FloatField(null=True)
    CB_CG_CD2 = models.FloatField(null=True)
    CB_CG_ND1 = models.FloatField(null=True)
    ND1_CG_CD2 = models.FloatField(null=True)
    CG_ND1_CE1 = models.FloatField(null=True)
    ND1_CE1_NE2 = models.FloatField(null=True)
    CE1_NE2_CD2 = models.FloatField(null=True)
    CG_CD2_NE2 = models.FloatField(null=True)


class Sidechain_ILE(models.Model):
    CB_CG1 = models.FloatField(null=True)
    CG1_CD1 = models.FloatField(null=True)
    CB_CG2 = models.FloatField(null=True)
    CA_CB_CG2 = models.FloatField(null=True)
    CA_CB_CG1 = models.FloatField(null=True)
    CG1_CB_CG2 = models.FloatField(null=True)
    CB_CG1_CD1 = models.FloatField(null=True)


class Sidechain_LEU(models.Model):
    CB_CG = models.FloatField(null=True)
    CG_CD1 = models.FloatField(null=True)
    CG_CD2 = models.FloatField(null=True)
    CA_CB_CG = models.FloatField(null=True)
    CB_CG_CD1 = models.FloatField(null=True)
    CB_CG_CD2 = models.FloatField(null=True)
    CD1_CG_CD2 = models.FloatField(null=True)


class Sidechain_LYS(models.Model):
    CB_CG = models.FloatField(null=True)
    CG_CD = models.FloatField(null=True)
    CD_CE = models.FloatField(null=True)
    CE_NZ = models.FloatField(null=True)
    CA_CB_CG = models.FloatField(null=True)
    CB_CG_CD = models.FloatField(null=True)
    CG_CD_CE = models.FloatField(null=True)
    CD_CE_NZ = models.FloatField(null=True)


class Sidechain_MET(models.Model):
    CB_CG = models.FloatField(null=True)
    CG_SD = models.FloatField(null=True)
    SD_CE = models.FloatField(null=True)
    CA_CB_CG = models.FloatField(null=True)
    CB_CG_SD = models.FloatField(null=True)
    CG_SD_CE = models.FloatField(null=True)


class Sidechain_PHE(models.Model):
    CB_CG = models.FloatField(null=True)
    CG_CD1 = models.FloatField(null=True)
    CG_CD2 = models.FloatField(null=True)
    CD1_CE1 = models.FloatField(null=True)
    CD2_CE2 = models.FloatField(null=True)
    CZ_CE1 = models.FloatField(null=True)
    CZ_CE2 = models.FloatField(null=True)
    CA_CB_CG = models.FloatField(null=True)
    CB_CG_CD1 = models.FloatField(null=True)
    CB_CG_CD2 = models.FloatField(null=True)
    CD1_CG_CD2 = models.FloatField(null=True)
    CE1_CZ_CE2 = models.FloatField(null=True)
    CG_CD1_CE1 = models.FloatField(null=True)
    CG_CD2_CE2 = models.FloatField(null=True)
    CZ_CE1_CD1 = models.FloatField(null=True)
    CZ_CE2_CD2 = models.FloatField(null=True)


class Sidechain_PRO(models.Model):
    CB_CG = models.FloatField(null=True)
    CG_CD = models.FloatField(null=True)
    CD_N = models.FloatField(null=True)
    CA_CB_CG = models.FloatField(null=True)
    CB_CG_CD = models.FloatField(null=True)
    CG_CD_N = models.FloatField(null=True)
    CD_N_CA = models.FloatField(null=True)
    CD_N_C_1 = models.FloatField(null=True)


class Sidechain_SER(models.Model):
    CB_OG = models.FloatField(null=True)
    CA_CB_OG = models.FloatField(null=True)


class Sidechain_THR(models.Model):
    CB_OG1 = models.FloatField(null=True)
    CB_CG2 = models.FloatField(null=True)
    CA_CB_OG1 = models.FloatField(null=True)
    CA_CB_CG2 = models.FloatField(null=True)
    OG1_CB_CG2 = models.FloatField(null=True)


class Sidechain_TRP(models.Model):
    CB_CG = models.FloatField(null=True)
    CG_CD1 = models.FloatField(null=True)
    CD1_NE1 = models.FloatField(null=True)
    NE1_CE2 = models.FloatField(null=True)
    CE2_CZ2 = models.FloatField(null=True)
    CZ2_CH2 = models.FloatField(null=True)
    CH2_CZ3 = models.FloatField(null=True)
    CZ3_CE3 = models.FloatField(null=True)
    CD2_CE3 = models.FloatField(null=True)
    CG_CD2 = models.FloatField(null=True)
    CD2_CE2 = models.FloatField(null=True)
    CA_CB_CG = models.FloatField(null=True)
    CB_CG_CD1 = models.FloatField(null=True)
    CB_CG_CD2 = models.FloatField(null=True)
    CD2_CG_CD1 = models.FloatField(null=True)
    CG_CD1_NE1 = models.FloatField(null=True)
    CD1_NE1_CE2 = models.FloatField(null=True)
    NE1_CE2_CD2 = models.FloatField(null=True)
    CG_CD2_CE2 = models.FloatField(null=True)
    CE2_CD2_CE3 = models.FloatField(null=True)
    CD2_CE2_CZ2 = models.FloatField(null=True)
    CE2_CZ2_CH2 = models.FloatField(null=True)
    CZ2_CH2_CZ3 = models.FloatField(null=True)
    CE3_CZ3_CH2 = models.FloatField(null=True)
    CD2_CE3_CZ3 = models.FloatField(null=True)
    CG_CD2_CE3 = models.FloatField(null=True)
    NE1_CE2_CZ2 = models.FloatField(null=True)


class Sidechain_TYR(models.Model):
    CB_CG = models.FloatField(null=True)
    CG_CD1 = models.FloatField(null=True)
    CG_CD2 = models.FloatField(null=True)
    CD1_CE1 = models.FloatField(null=True)
    CD2_CE2 = models.FloatField(null=True)
    CZ_CE1 = models.FloatField(null=True)
    CZ_CE2 = models.FloatField(null=True)
    CZ_OH = models.FloatField(null=True)
    CA_CB_CG = models.FloatField(null=True)
    CB_CG_CD1 = models.FloatField(null=True)
    CB_CG_CD2 = models.FloatField(null=True)
    CD1_CG_CD2 = models.FloatField(null=True)
    CE1_CZ_CE2 = models.FloatField(null=True)
    CG_CD1_CE1 = models.FloatField(null=True)
    CG_CD2_CE2 = models.FloatField(null=True)
    CZ_CE1_CD1 = models.FloatField(null=True)
    CZ_CE2_CD2 = models.FloatField(null=True)
    CE1_CZ_OH = models.FloatField(null=True)
    CE2_CZ_OH = models.FloatField(null=True)


class Sidechain_VAL(models.Model):
    CB_CG1 = models.FloatField(null=True)
    CB_CG2 = models.FloatField(null=True)
    CA_CB_CG1 = models.FloatField(null=True)
    CA_CB_CG2 = models.FloatField(null=True)
    CG1_CB_CG2 = models.FloatField(null=True)


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
    prev            = models.ForeignKey('self', related_name='prev_next', null=True)
    next            = models.ForeignKey('self', related_name='next_prev', null=True)
    aa              = models.CharField(max_length=1, choices=AA_CHOICES) # new type
    chainID         = models.CharField(max_length=1) # integer id
    oldID           = models.CharField(max_length=5, null=True)# id[icode] from pdb file
    chainIndex      = models.PositiveIntegerField()
    a1              = models.FloatField(null=True)
    a2              = models.FloatField(null=True)
    a3              = models.FloatField()
    a4              = models.FloatField(null=True)
    a5              = models.FloatField()
    a6              = models.FloatField(null=True)
    a7              = models.FloatField(null=True)
    L1              = models.FloatField(null=True)
    L2              = models.FloatField()
    L3              = models.FloatField(null=True)
    L4              = models.FloatField()
    L5              = models.FloatField()
    ss              = models.CharField(max_length=1, choices=SS_CHOICES) # new type (was blob, but all entries 1 char)
    phi             = models.FloatField(null=True)
    psi             = models.FloatField(null=True)
    ome             = models.FloatField(null=True)
    omep           = models.FloatField(null=True)
    chi1            = models.FloatField(null=True)
    chi2            = models.FloatField(null=True)
    chi3            = models.FloatField(null=True)
    chi4            = models.FloatField(null=True)
    chi5            = models.FloatField(null=True)
    bm              = models.FloatField()
    bs              = models.FloatField()
    bg              = models.FloatField(null=True)
    h_bond_energy   = models.FloatField()
    zeta            = models.FloatField(null=True)
    terminal_flag   = models.BooleanField(default=False)#indicates this residue is next to a chain break
    xpr             = models.BooleanField(default=False) # this field may not be necessary; it has never been implemented

    sidechain_ARG = models.OneToOneField(Sidechain_ARG, related_name="residue", null=True)
    sidechain_ASN = models.OneToOneField(Sidechain_ASN, related_name="residue", null=True)
    sidechain_ASP = models.OneToOneField(Sidechain_ASP, related_name="residue", null=True)
    sidechain_CYS = models.OneToOneField(Sidechain_CYS, related_name="residue", null=True)
    sidechain_GLN = models.OneToOneField(Sidechain_GLN, related_name="residue", null=True)
    sidechain_GLU = models.OneToOneField(Sidechain_GLU, related_name="residue", null=True)
    sidechain_HIS = models.OneToOneField(Sidechain_HIS, related_name="residue", null=True)
    sidechain_ILE = models.OneToOneField(Sidechain_ILE, related_name="residue", null=True)
    sidechain_LEU = models.OneToOneField(Sidechain_LEU, related_name="residue", null=True)
    sidechain_LYS = models.OneToOneField(Sidechain_LYS, related_name="residue", null=True)
    sidechain_MET = models.OneToOneField(Sidechain_MET, related_name="residue", null=True)
    sidechain_PHE = models.OneToOneField(Sidechain_PHE, related_name="residue", null=True)
    sidechain_PRO = models.OneToOneField(Sidechain_PRO, related_name="residue", null=True)
    sidechain_SER = models.OneToOneField(Sidechain_SER, related_name="residue", null=True)
    sidechain_THR = models.OneToOneField(Sidechain_THR, related_name="residue", null=True)
    sidechain_TRP = models.OneToOneField(Sidechain_TRP, related_name="residue", null=True)
    sidechain_TYR = models.OneToOneField(Sidechain_TYR, related_name="residue", null=True)
    sidechain_VAL = models.OneToOneField(Sidechain_VAL, related_name="residue", null=True)

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
