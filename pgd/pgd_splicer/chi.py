"""
Dictionary used to map peptides to the list of atoms used to make up its chi's
"""

PROTEIN_ORDER = ['VAL','LEU','ILE','PRO','PHE','TYR','TRP','MET','CYS','SER','THR','ASP','GLU','HIS','LYS','ARG','ASN','GLN']

CHI_MAP = {
    'VAL':(('N', 'CA', 'CB', 'CG1'),),
    'LEU':(
            ('N',   'CA',  'CB',  'CG'),
            ('CA',  'CB',  'CG',  'CD1'),
        ),
    'ILE':(
        ('N',   'CA',  'CB',   'CG1'),
        ('CA',  'CB',  'CG1',  'CD1')
        ),
    'PRO':(('N',  'CA', 'CB', 'CG'),
           ('CA', 'CB', 'CG', 'CD'),
           ('CB', 'CG', 'CD', 'N'),
           ('CG', 'CD', 'N',  'CA')
        ),
    'PHE':(
        ('N',   'CA',  'CB',  'CG'),
        ('CA',  'CB',  'CG',  'CD1')
        ),
    'TYR':(
        ('N',   'CA',  'CB',  'CG'),
        ('CA',  'CB',  'CG',  'CD1')
        ),
    'TRP':(
        ('N',   'CA',  'CB',  'CG'),
        ('CA',  'CB',  'CG',  'CD1')
        ),
    'MET':(
        ('N',   'CA',  'CB',  'CG'),
        ('CA',  'CB',  'CG',  'SD'),
        ('CB',  'CG',  'SD',  'CE')
        ),
    'CYS':(('N', 'CA', 'CB', 'SG'),),
    'SER':(('N', 'CA', 'CB', 'OG'),),
    'THR':(('N', 'CA', 'CB', 'OG1'),),
    'ASP':(
        ('N',   'CA',  'CB',  'CG'),
        ('CA',  'CB',  'CG',  'OD1')
        ),
    'GLU':(
        ('N',   'CA',  'CB',  'CG'),
        ('CA',  'CB',  'CG',  'CD'),
        ('CB',  'CG',  'CD',  'OE1')
        ),
    'HIS':(
        ('N',   'CA',  'CB',  'CG'),
        ('CA',  'CB',  'CG',  'ND'),
        ('CA','CB','CD','ND1')
        ),
    'LYS':(
        ('N',   'CA',  'CB',  'CG'),
        ('CA',  'CB',  'CG',  'CD'),
        ('CB',  'CG',  'CD',  'CE'),
        ('CG',  'CD',  'CE',  'NZ')
        ),
    'ARG':(
        ('N',   'CA',  'CB',  'CG'),
        ('CA',  'CB',  'CG',  'CD'),
        ('CB',  'CG',  'CD',  'NE'),
        ('CG',  'CD',  'NE',  'CZ'),
        ('CD',  'NE',  'CZ',  'NH1')
        ),
    'ASN':(
        ('N',   'CA',  'CB',  'CG'),
        ('CA',  'CB',  'CG',  'OD1')
        ),
    'GLN':(
        ('N',   'CA',  'CB',  'CG'),
        ('CA',  'CB',  'CG',  'CD'),
        ('CB',  'CG',  'CD',  'OE1')
        )
    }


# CHI corrections are symettrical values that aren't guarunteed to be recorded
# in the correct order by within the PDB file.  The correct order can be
# determined by checking for the larger of two Chi values.  If the 2nd value
# listed in CHI_CORRECTIONS_TESTS is larger, then corrections are needed. There
# may be multiple corrections per Residue which are listed in CHI_CORRECTIONS.
# each pair of atoms will be swapped so that any future calculation will use
# the correct values.
#
# See: Ticket #1545
CHI_CORRECTIONS_TESTS = {
    'ASP':[('CA',  'CB',  'CG',  'OD1'), ('CA',  'CB',  'CG',  'OD2')],
    'GLU':[('CB',  'CG',  'CD',  'OE1'), ('CB',  'CG',  'CD',  'OE2')],
    'ARG':[('CD',  'NE',  'CZ',  'NH1'), ('CD',  'NE',  'CZ',  'NH2')],
    'PHE':[('CA',  'CB',  'CG',  'CD1'), ('CA',  'CB',  'CG',  'CD2')],
    'TYR':[('CA',  'CB',  'CG',  'CD1'), ('CA',  'CB',  'CG',  'CD2')],
}

CHI_CORRECTIONS = {
    'ASP':[('OD1','OD2')],
    'GLU':[('OE1','OE2')],
    'ARG':[('NH1','NH2')],
    'PHE':[('CD1','CD2'), ('CE1','CE2')],
    'TYR':[('CD1','CD2'), ('CE1','CE2')],
}