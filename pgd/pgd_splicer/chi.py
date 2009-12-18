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
    'PRO':(('N', 'CA', 'CB', 'CG'),),
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
        ('CA',  'CB',  'CG',  'OE1')
        ),
    'HIS':(
        ('N',   'CA',  'CB',  'CG'),
        ('CA',  'CB',  'CG',  'ND')
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
        ('CG',  'CD',  'NE',  'CZ')
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
