from django.db import models
from django.contrib.auth.models import User
from pgd_constants import AA_CHOICES, SS_CHOICES
from SearchParser import parse_search

# Search
# A search query submitted by a user
class Search(models.Model):
    user             = models.ForeignKey(User, null=True)
    codes_include    = models.BooleanField(null=True)
    threshold        = models.IntegerField(null=True)
    resolution_min   = models.FloatField(null=True)
    resolution_max   = models.FloatField(null=True)
    _querySet = None

    # returns the query set that represents this search
    def querySet(self):
        #create querySet if not needed
        if not self._querySet:
            self._querySet = parse_search(self)

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


    def __init__(self, *args, **kwargs):
        models.Model.__init__(self, *args, **kwargs)

        # populate 'aa' with a dictionary of allowed values from AA_CHOICES
        self.aa = dict([(j[1],1 if self.aa_int == None else 1&self.aa_int>>i) for i,j in enumerate(AA_CHOICES)])