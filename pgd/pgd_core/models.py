from django.db import models
from constants import AA_CHOICES, SS_CHOICES, Subscripter

# Protein model
# (was 'protein_info')
# Contains the general information about a protein
class Protein(models.Model):
    code        = models.CharField(max_length=4, primary_key=True, unique=True)
    threshold   = models.IntegerField() # new type; this should probably be a boolean type
    resolution  = models.FloatField()
    rfactor     = models.FloatField()

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

# Residue model
# (was 'protein')
# Contains information regarding each amino acid and its geometry
# (Note: fields need to be commented)
class Residue(models.Model):
    chain           = models.ForeignKey(Chain, related_name='residues')
    protein         = models.ForeignKey(Protein, related_name='residues')
    aa              = models.CharField(max_length=1, choices=AA_CHOICES) # new type
    chainID         = models.CharField(max_length=1)
    chainIndex      = models.PositiveIntegerField()
    a1              = models.FloatField()
    a2              = models.FloatField()
    a3              = models.FloatField()
    a4              = models.FloatField()
    a5              = models.FloatField()
    a6              = models.FloatField()
    a7              = models.FloatField()
    L1              = models.FloatField()
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
    bg              = models.FloatField()
    h_bond_energy   = models.FloatField()
    zeta            = models.FloatField()
    terminal_flag   = models.BooleanField() 
    xpr             = models.BooleanField() # this field may not be necessary; it has never been implemented

    def __str__(self):
        return '%d' % self.newID