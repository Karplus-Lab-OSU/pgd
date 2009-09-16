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
        ['CD', 'N', 'CA']],
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
