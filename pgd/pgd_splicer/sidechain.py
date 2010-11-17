bond_lengths = {
    'CYS': [
        ['CB', 'SG']
        ],
    'GLN': [
        ['CB', 'CG'],
        ['CG', 'CD'],
        ['CD', 'OE1'],
        ['CD', 'NE2']
        ],
    'ILE': [
        ['CB', 'CG1'],
        ['CG1', 'CD1'],
        ['CB', 'CG2']],
    'SER': [
        ['CB', 'OG']],
    'VAL': [
        ['CB', 'CG1'],
        ['CB', 'CG2']],
    'LYS': [
        ['CB', 'CG'],
        ['CG', 'CD'],
        ['CD', 'CE'],
        ['CE', 'NZ']],
    'PRO': [
        ['CB', 'CG'],
        ['CG', 'CD'],
        ['CD', 'N']],
    'THR': [
        ['CB', 'OG1'],
        ['CB', 'CG2']],
    'PHE': [
        ['CB', 'CG'],
        ['CG', 'CD1'],
        ['CG', 'CD2'],
        ['CD1', 'CE1'],
        ['CD2', 'CE2'],
        ['CZ', 'CE1'],
        ['CZ', 'CE2']],
    'GLU': [
        ['CB', 'CG'],
        ['CG', 'CD'],
        ['CD', 'OE1'],
        ['CD', 'OE2']],
    'HIS': [
        ['CB', 'CG'],
        ['CG', 'ND1'],
        ['ND1', 'CE1'],
        ['CE1', 'NE2'],
        ['NE2', 'CD2'],
        ['CG', 'CD2']],
    'MET': [
        ['CB', 'CG'],
        ['CG', 'SD'],
        ['SD', 'CE']],
    'ASP': [
        ['CB', 'CG'],
        ['CG', 'OD1'],
        ['CG', 'OD2']],
    'LEU': [
        ['CB', 'CG'],
        ['CG', 'CD1'],
        ['CG', 'CD2']],
    'ARG': [
        ['CB', 'CG'],
        ['CG', 'CD'],
        ['CD', 'NE'],
        ['NE', 'CZ'],
        ['CZ', 'NH1'],
        ['CZ', 'NH2']],
    'TRP': [
        ['CB', 'CG'],
        ['CG', 'CD1'],
        ['CD1', 'NE1'],
        ['NE1', 'CE2'],
        ['CE2', 'CZ2'],
        ['CZ2', 'CH2'],
        ['CH2', 'CZ3'],
        ['CZ3', 'CE3'],
        ['CD2', 'CE3'],
        ['CG', 'CD2'],
        ['CD2', 'CE2']],
    'ASN': [
        ['CB', 'CG'],
        ['CG', 'OD1'],
        ['CG', 'ND2']],
    'TYR': [
        ['CB', 'CG'],
        ['CG', 'CD1'],
        ['CG', 'CD2'],
        ['CD1', 'CE1'],
        ['CD2', 'CE2'],
        ['CZ', 'CE1'],
        ['CZ', 'CE2'],
        ['CZ', 'OH']]
    }


bond_angles = {
    'CYS': [
        ['CA', 'CB', 'SG']],
    'GLN': [
        ['CA', 'CB', 'CG'],
        ['CB', 'CG', 'CD'],
        ['CG', 'CD', 'OE1'],
        ['CG', 'CD', 'NE2'],
        ['OE1', 'CD', 'NE2']],
    'ILE': [
        ['CA', 'CB', 'CG2'],
        ['CA', 'CB', 'CG1'],
        ['CG1', 'CB', 'CG2'],
        ['CB', 'CG1', 'CD1']],
    'SER': [
        ['CA', 'CB', 'OG']],
    'VAL': [
        ['CA', 'CB', 'CG1'],
        ['CA', 'CB', 'CG2'],
        ['CG1', 'CB', 'CG2']],
    'LYS': [
        ['CA', 'CB', 'CG'],
        ['CB', 'CG', 'CD'],
        ['CG', 'CD', 'CE'],
        ['CD', 'CE', 'NZ']],
    'PRO': [
        ['CA', 'CB', 'CG'],
        ['CB', 'CG', 'CD'],
        ['CG', 'CD', 'N'],
        ['CD', 'N', 'CA'],
        ['CD', 'N', 'C-1']],
    'THR': [
        ['CA', 'CB', 'OG1'],
        ['CA', 'CB', 'CG2'],
        ['OG1', 'CB', 'CG2']], 
    'PHE': [
        ['CA', 'CB', 'CG'],
        ['CB', 'CG', 'CD1'],
        ['CB', 'CG', 'CD2'],
        ['CD1', 'CG', 'CD2'],
        ['CE1', 'CZ', 'CE2'],
        ['CG', 'CD1', 'CE1'],
        ['CG', 'CD2', 'CE2'],
        ['CZ', 'CE1', 'CD1'],
        ['CZ', 'CE2', 'CD2']],
    'GLU': [
        ['CA', 'CB', 'CG'],
        ['CB', 'CG', 'CD'],
        ['CG', 'CD', 'OE1'],
        ['CG', 'CD', 'OE2'],
        ['OE1', 'CD', 'OE2']],
    'HIS': [
        ['CA', 'CB', 'CG'],
        ['CB', 'CG', 'CD2'],
        ['CB', 'CG', 'ND1'],
        ['ND1', 'CG', 'CD2'],
        ['CG', 'ND1', 'CE1'],
        ['ND1', 'CE1', 'NE2'],
        ['CE1', 'NE2', 'CD2'],
        ['CG', 'CD2', 'NE2']],
    'MET': [
        ['CA', 'CB', 'CG'],
        ['CB', 'CG', 'SD'],
        ['CG', 'SD', 'CE']],
    'ASP': [
        ['CA', 'CB', 'CG'],
        ['CB', 'CG', 'OD1'],
        ['CB', 'CG', 'OD2'],
        ['OD1', 'CG', 'OD2']],
    'LEU': [
        ['CA', 'CB', 'CG'],
        ['CB', 'CG', 'CD1'],
        ['CB', 'CG', 'CD2'],
        ['CD1', 'CG', 'CD2']],
    'ARG': [
        ['CA', 'CB', 'CG'],
        ['CB', 'CG', 'CD'],
        ['CG', 'CD', 'NE'],
        ['CD', 'NE', 'CZ'],
        ['NE', 'CZ', 'NH1'],
        ['NE', 'CZ', 'NH2'],
        ['NH1', 'CZ', 'NH2']],
    'TRP': [
        ['CA', 'CB', 'CG'],
        ['CB', 'CG', 'CD1'],
        ['CB', 'CG', 'CD2'],
        ['CD2', 'CG', 'CD1'],
        ['CG', 'CD1', 'NE1'],
        ['CD1', 'NE1', 'CE2'],
        ['NE1', 'CE2', 'CD2'],
        ['CG', 'CD2', 'CE2'],
        ['CE2', 'CD2', 'CE3'],
        ['CD2', 'CE2', 'CZ2'],
        ['CE2', 'CZ2', 'CH2'],
        ['CZ2', 'CH2', 'CZ3'],
        ['CE3', 'CZ3', 'CH2'],
        ['CD2', 'CE3', 'CZ3'],
        ['CG', 'CD2', 'CE3'],
        ['NE1', 'CE2', 'CZ2']],
    'ASN': [
        ['CA', 'CB', 'CG'],
        ['CB', 'CG', 'OD1'],
        ['CB', 'CG', 'ND2'],
        ['OD1', 'CG', 'ND2']],
    'TYR': [
        ['CA', 'CB', 'CG'],
        ['CB', 'CG', 'CD1'],
        ['CB', 'CG', 'CD2'],
        ['CD1', 'CG', 'CD2'],
        ['CE1', 'CZ', 'CE2'],
        ['CG', 'CD1', 'CE1'],
        ['CG', 'CD2', 'CE2'],
        ['CZ', 'CE1', 'CD1'],
        ['CZ', 'CE2', 'CD2'],
        ['CE1', 'CZ', 'OH'],
        ['CE2', 'CZ', 'OH']]
    }

bond_angles_string_dict = {
    'CYS': [
        'CA_CB_SG'],
    'GLN': [
        'CA_CB_CG',
        'CB_CG_CD',
        'CG_CD_OE1',
        'CG_CD_NE2',
        'OE1_CD_NE2'],
    'ILE': [
        'CA_CB_CG2',
        'CA_CB_CG1',
        'CG1_CB_CG2',
        'CB_CG1_CD1'],
    'SER': [
        'CA_CB_OG'],
    'VAL': [
        'CA_CB_CG1',
        'CA_CB_CG2',
        'CG1_CB_CG2'],
    'LYS': [
        'CA_CB_CG',
        'CB_CG_CD',
        'CG_CD_CE',
        'CD_CE_NZ'],
    'PRO': [
        'CA_CB_CG',
        'CB_CG_CD',
        'CG_CD_N',
        'CD_N_CA',
        'CD_N_C-1'],
    'THR': [
        'CA_CB_OG1',
        'CA_CB_CG2',
        'OG1_CB_CG2'], 
    'PHE': [
        'CA_CB_CG',
        'CB_CG_CD1',
        'CB_CG_CD2',
        'CD1_CG_CD2',
        'CE1_CZ_CE2',
        'CG_CD1_CE1',
        'CG_CD2_CE2',
        'CZ_CE1_CD1',
        'CZ_CE2_CD2'],
    'GLU': [
        'CA_CB_CG',
        'CB_CG_CD',
        'CG_CD_OE1',
        'CG_CD_OE2',
        'OE1_CD_OE2'],
    'HIS': [
        'CA_CB_CG',
        'CB_CG_CD2',
        'CB_CG_ND1',
        'ND1_CG_CD2',
        'CG_ND1_CE1',
        'ND1_CE1_NE2',
        'CE1_NE2_CD2',
        'CG_CD2_NE2'],
    'MET': [
        'CA_CB_CG',
        'CB_CG_SD',
        'CG_SD_CE'],
    'ASP': [
        'CA_CB_CG',
        'CB_CG_OD1',
        'CB_CG_OD2',
        'OD1_CG_OD2'],
    'LEU': [
        'CA_CB_CG',
        'CB_CG_CD1',
        'CB_CG_CD2',
        'CD1_CG_CD2'],
    'ARG': [
        'CA_CB_CG',
        'CB_CG_CD',
        'CG_CD_NE',
        'CD_NE_CZ',
        'NE_CZ_NH1',
        'NE_CZ_NH2',
        'NH1_CZ_NH2'],
    'TRP': [
        'CA_CB_CG',
        'CB_CG_CD1',
        'CB_CG_CD2',
        'CD2_CG_CD1',
        'CG_CD1_NE1',
        'CD1_NE1_CE2',
        'NE1_CE2_CD2',
        'CG_CD2_CE2',
        'CE2_CD2_CE3',
        'CD2_CE2_CZ2',
        'CE2_CZ2_CH2',
        'CZ2_CH2_CZ3',
        'CE3_CZ3_CH2',
        'CD2_CE3_CZ3',
        'CG_CD2_CE3',
        'NE1_CE2_CZ2'],
    'ASN': [
        'CA_CB_CG',
        'CB_CG_OD1',
        'CB_CG_ND2',
        'OD1_CG_ND2'],
    'TYR': [
        'CA_CB_CG',
        'CB_CG_CD1',
        'CB_CG_CD2',
        'CD1_CG_CD2',
        'CE1_CZ_CE2',
        'CG_CD1_CE1',
        'CG_CD2_CE2',
        'CZ_CE1_CD1',
        'CZ_CE2_CD2',
        'CE1_CZ_OH',
        'CE2_CZ_OH']
    }
bond_lengths_string_dict = {
    'CYS': [
        'CB_SG'
        ],
    'GLN': [
        'CB_CG',
        'CG_CD',
        'CD_OE1',
        'CD_NE2'
        ],
    'ILE': [
        'CB_CG1',
        'CG1_CD1',
        'CB_CG2'],
    'SER': [
        'CB_OG'],
    'VAL': [
        'CB_CG1',
        'CB_CG2'],
    'LYS': [
        'CB_CG',
        'CG_CD',
        'CD_CE',
        'CE_NZ'],
    'PRO': [
        'CB_CG',
        'CG_CD',
        'CD_N'],
    'THR': [
        'CB_OG1',
        'CB_CG2'],
    'PHE': [
        'CB_CG',
        'CG_CD1',
        'CG_CD2',
        'CD1_CE1',
        'CD2_CE2',
        'CZ_CE1',
        'CZ_CE2'],
    'GLU': [
        'CB_CG',
        'CG_CD',
        'CD_OE1',
        'CD_OE2'],
    'HIS': [
        'CB_CG',
        'CG_ND1',
        'ND1_CE1',
        'CE1_NE2',
        'NE2_CD2',
        'CG_CD2'],
    'MET': [
        'CB_CG',
        'CG_SD',
        'SD_CE'],
    'ASP': [
        'CB_CG',
        'CG_OD1',
        'CG_OD2'],
    'LEU': [
        'CB_CG',
        'CG_CD1',
        'CG_CD2'],
    'ARG': [
        'CB_CG',
        'CD_NE',
        'CG_CD',
        'NE_CZ',
        'CZ_NH1',
        'CZ_NH2'],
    'TRP': [
        'CB_CG',
        'CG_CD1',
        'CD1_NE1',
        'NE1_CE2',
        'CE2_CZ2',
        'CZ2_CH2',
        'CH2_CZ3',
        'CZ3_CE3',
        'CD2_CE3',
        'CG_CD2',
        'CD2_CE2'],
    'ASN': [
        'CB_CG',
        'CG_OD1',
        'CG_ND2'],
    'TYR': [
        'CB_CG',
        'CG_CD1',
        'CG_CD2',
        'CD1_CE1',
        'CD2_CE2',
        'CZ_CE1',
        'CZ_CE2',
        'CZ_OH']
    }

# combined dictionary of angles and lengths, grouped by aa
sidechain_string_dict = {}
sidechain_string_dict.update(bond_angles_string_dict)
for k,v in bond_lengths_string_dict.items():
    sidechain_string_dict[k] = sidechain_string_dict[k] + v

# list of sidechain properties without aa type
sidechain_angle_list = ['%s_%s'%(k,v) for k, a in bond_angles_string_dict.items() for v in a ]
sidechain_length_list = ['%s_%s'%(k,v) for k, a in bond_lengths_string_dict.items() for v in a ]

# list of sidechain properties with aa type
sidechain_angle_relationship_list = ['%s__%s'%(k,v) for k, a in bond_angles_string_dict.items() for v in a ]
sidechain_length_relationship_list = ['%s__%s'%(k,v) for k, a in bond_lengths_string_dict.items() for v in a ]


aa_list = [
    'ARG',
    'ASN',
    'ASP',
    'CYS',
    'GLN',
    'GLU',
    'HIS',
    'ILE',
    'LEU',
    'LYS',
    'MET',
    'PHE',
    'PRO',
    'SER',
    'THR',
    'TRP',
    'TYR',
    'VAL'
    ]

if __name__ == '__main__':
    pass
    
