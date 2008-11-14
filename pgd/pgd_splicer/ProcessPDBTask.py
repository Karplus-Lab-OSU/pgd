#!/usr/bin/env python

import math, sys, mmLib.Structure, mmLib.FileIO, mmLib.AtomMath
import Bio.PDB
import shutil
import os

# Structure.AminoAcidResidue:
#
# calc_mainchain_bond_angle(self)
# Calculates main chain bond angles (N-CA-C, N-CA-CB, CB-CA-C, CA-C-O,
# CA-C-(next)N, C-(next residue)N-(next residue)CA) and returns the result as
# a 6-tuple in that order.
#
# calc_mainchain_bond_length(self)
# Calculates the main chain bond lengths: (N-CA, CA-C, C-O, CA-CB, CA-(next)N).
#
# calc_torsion_phi(self)
# Calculates the Phi torsion angle of the amino acid.
#
# calc_torsion_psi(self)
# Calculates the Psi torsion angle of the amino acid.
#
# calc_torsion_omega(self)
# Calculates the Omega torsion angle of the amino acid.

NO_VALUE = 999.9

def processPDB(file):
    bioPythonProps = parseWithBioPython(file)
    #print 'props: %s' % bioPythonProps
    return parseWithMMLib(file, bioPythonProps)


def decompressFile(src):
    tempfile = 'tmp/%s' % src
    dest = 'tmp/%s' % src[:-2]

    try:
        # copy the file to the tmp directory
        shutil.copyfile(src, tempfile)

        # decompress using unix decompress.  The files are stored used LZ.
        # algorithms i've tried will not decompress the files properly.  for
        # now use the linux decompress
        os.system('uncompress %s' % tempfile)

    except:
        #clean up resulting file on errors
        if os.path.exists(dest):
            os.path.remove(dest)

        return False

    finally:
        # clean up temp file no matter what
        if os.path.exists(tempfile):
            os.path.remove(tempfile)

    return dest

def parseWithBioPython(file):
    residueProps = {}
    decompressedFile = None
    tmp = './tmp'

    try:
        #create tmp workspace
        if os.path.exists(tmp):
            ownTempDir = False
        else:
            ownTempDir = True
            os.mkdir(tmp)

        #prep and open file
        decompressedFile = decompressFile(file)

        if decompressedFile:

            structure = Bio.PDB.PDBParser().get_structure('pdbname', decompressedFile)

            # dssp can't do multiple models. if we ever need to, we'll have to 
            # iterate through them
            dssp = Bio.PDB.DSSP(model=structure[0], pdb_file=decompressedFile, dssp='dsspcmbi')

            #iterate residues
            for res in structure[0].get_residues():
                hetflag, resseq, icode = res.get_id()
                chain = res.get_parent().get_id()
                if hetflag != ' ':
                    continue

                # get properties using dssp
                residue, secondary_structure, accessibility, relative_accessibility = dssp[(chain, resseq)]

                # save in dictionary for later use
                residueProps[resseq] = {'ss':secondary_structure}
    finally:
        #clean up any files in tmp directory no matter what
        if decompressedFile and os.path.exists(decompressedFile):
            os.remove(decompressedFile)

        if ownTempDir and os.path.exists(tmp):
            os.removedirs(tmp)

    return residueProps

def parseWithMMLib(file, props):
    struct = mmLib.FileIO.LoadStructure(file=file)
    lowerCaseIndicators = ['H','G','E','T']
    ret = {}
    for r in struct.iter_amino_acids():
        # get offsets for residues left and right of this one
        next_res = r.get_offset_residue(1)
        prev_res = r.get_offset_residue(-1)

        #Main Bond Angle
        r.a3, r.a2, r.a4, r.a5, a6, a1 = \
              r.calc_mainchain_bond_angle()

        #Main Bond Length
        r.l2, r.l4, r.l5, r.l3, r.l6 = \
              r.calc_mainchain_bond_length()
        r.l7 = NO_VALUE


        # if there is a residue before this one, get the angles
        if prev_res:
            r.l1 = prev_res.l6
        else:
            # n-terminus
            r.l1 = NO_VALUE
            r.a1 = math.radians(NO_VALUE)
            r.a6 = math.radians(NO_VALUE)
            r.a7 = math.radians(NO_VALUE)


        #if there is a residue after this one, get properties for that relationship
        if next_res:
            next_res.a1 = a1
            next_res.a6 = a6
            aN = r.get_atom('N')
            aCA = r.get_atom('CA')
            aO = r.get_atom('O')
            aC = r.get_atom('C')
            naN = next_res.get_atom('N')

            next_res.l2 = mmLib.AtomMath.calc_distance(aN, aCA)
            r.l7 = next_res.l2
            next_res.a7 = mmLib.AtomMath.calc_angle(aO, aC, naN)
        else:
            # c-terminus
            r.l7 = NO_VALUE

        r.x1, r.x2, r.x3, r.x4 = r.calc_torsion_chi()
        r.phi = r.calc_torsion_phi()
        r.psi = r.calc_torsion_psi()
        r.ome = r.calc_torsion_omega()

        if not r.l3:
            # glycine
            r.l3 = NO_VALUE
            r.a2 = math.radians(NO_VALUE)
            r.a4 = math.radians(NO_VALUE)
        if not r.l6:
            # c-terminus
            r.l6 = NO_VALUE
        if not r.phi:
            r.phi = math.radians(NO_VALUE)
        if not r.psi:
            r.psi = math.radians(NO_VALUE)
        if not r.ome:
            r.ome = math.radians(NO_VALUE)

        # check for ends of consecutive groupings.  ends must be lowercase if they are in [hget]
        ss = props[int(r.fragment_id)]['ss']
        if ((prev_res == None or next_res == None)) and ss in lowerCaseIndicators:
            ss = ss.lower()
            print 'lower!'
        else:
            print '%s : %s %s %s' % (ss, prev_res, next_res, ss in lowerCaseIndicators)


        # This accounts for the possibility of missing atoms by initializing
        # all of the values to NO_VALUE
        dihedral_list = ['phi', 'psi', 'ome']
        angle_list = ['a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7']
        length_list = ['l1', 'l2', 'l3', 'l4', 'l5', 'l6', 'l7']
        r = initialize_geometry(r, dihedral_list, 'angle')
        r = initialize_geometry(r, angle_list, 'angle')
        r = initialize_geometry(r, length_list, 'length')

        if not r.x1:
            r.x1 = math.radians(NO_VALUE)
        if not r.x2:
            r.x2 = math.radians(NO_VALUE)
        if not r.x3:
            r.x3 = math.radians(NO_VALUE)
        if not r.x4:
            r.x4 = math.radians(NO_VALUE)

        if not r.phi:
            r.phi = math.radians(NO_VALUE)
        if not r.psi:
            r.psi = math.radians(NO_VALUE)
        if not r.ome:
            r.ome = math.radians(NO_VALUE)

        r.zeta = mmLib.AtomMath.calc_torsion_angle(
                   r.get_atom('CA'), 
                   r.get_atom('N'), 
                   r.get_atom('C'), 
                   r.get_atom('CB')
                   )
        if r.zeta is None:
            r.zeta = math.radians(NO_VALUE)

        all_atoms = mmLib.Structure.AtomList()
        mainchain_atoms = mmLib.Structure.AtomList()
        sidechain_atoms = mmLib.Structure.AtomList()

        for atom in r:
            if atom.name in ('N', 'CA', 'C', 'O'):
                mainchain_atoms.append(atom)
            else:
                # ignore hydrogens
                if atom.name == 'H':
                    continue
                sidechain_atoms.append(atom)
                if atom.name == 'CG':
                    try:
                        r.b_gamma
                        print "overwriting b_gamma for", r
                    except:
                        pass
                    r.b_gamma = atom.temp_factor

        if mainchain_atoms:
            r.b_main = mainchain_atoms.calc_adv_temp_factor()
        else:
            r.b_main = NO_VALUE
        if sidechain_atoms:
            r.b_side = sidechain_atoms.calc_adv_temp_factor()
        else:
            r.b_side = NO_VALUE
        try:
            r.b_gamma
        except:
            r.b_gamma = NO_VALUE

        
        


        r.props = {
            'name': r.res_name,
            'id': r.fragment_id,
            'phi': math.degrees(r.phi),
            'psi': math.degrees(r.psi),
            'ome': math.degrees(r.ome),
            'l1': r.l1,
            'l2': r.l2,
            'l3': r.l3,
            'l4': r.l4,
            'l5': r.l5,
            'l6': r.l6,
            'l7': r.l7,
            'a1': math.degrees(r.a1),
            'a2': math.degrees(r.a2),
            'a3': math.degrees(r.a3),
            'a4': math.degrees(r.a4),
            'a5': math.degrees(r.a5),
            'a6': math.degrees(r.a6),
            'a7': math.degrees(r.a7),
            'ss': ss,
            'chain': r.chain_id,
            'hb': 0.00,
            'zeta': math.degrees(r.zeta),
            'bg': r.b_gamma,
            'bm': r.b_main,
            'bs': r.b_side,
            'zero': 0.00,
            'x1': math.degrees(r.x1),
            'x2': math.degrees(r.x2),
            'x3': math.degrees(r.x3),
            'x4': math.degrees(r.x4),
            }
#        print '%(name)s%(id)6s%(phi)7.1f%(psi)6.1f%(ome)6.1f' % r.props,
        #print '%(name)s%(id)6s%(a6)7.1f%(a7)6.1f%(a1)6.1f' % r.props,
        #print '%(l1)5.3f%(l2)6.3f%(l4)6.3f%(l5)6.3f%(a3)6.1f%(a5)6.1f%(a2)6.1f%(a4)6.1f%(l3)6.3f' % r.props,
        #print
        #print r.props
        ret[int(r.fragment_id)]= r.props
    return ret

def initialize_geometry(residue, geometry_list, type):
    for item in geometry_list:
        if (getattr(residue, item)) is None:
            if type == 'angle':
                setattr(residue, item, math.radians(NO_VALUE))
            elif type == 'length':
                setattr(residue, item, NO_VALUE)
            else:
                print "Don't know how to deal with type", type
    return residue

if __name__ == '__main__':
    props = processPDB(sys.argv[1])
    s = ''
    for key in props:
        s = '%s%s' % (s, props[key]['ss'])
    print s

