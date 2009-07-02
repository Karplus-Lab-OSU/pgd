#!/usr/bin/env python
if __name__ == '__main__':
    import sys
    import os

    #python magic to add the current directory to the pythonpath
    sys.path.append(os.getcwd())

    # ==========================================================
    # Setup django environment 
    # ==========================================================
    if not os.environ.has_key('DJANGO_SETTINGS_MODULE'):
        os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'
    # ==========================================================
    # Done setting up django environment
    # ==========================================================

from tasks.tasks import *

from pgd_splicer.models import *
from pgd_core.models import Protein as ProteinModel
from pgd_core.models import Chain as ChainModel
from pgd_core.models import Residue as ResidueModel

from django.db import transaction

import math, sys, mmLib.Structure, mmLib.FileIO, mmLib.AtomMath
import Bio.PDB
from Bio.PDB import *
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

AA3to1 =  {
    'ALA' : 'a',
    'ARG' : 'r',
    'ASN' : 'n',
    'ASP' : 'd',
    'CYS' : 'c',
    'GLU' : 'e',
    'GLN' : 'q',
    'GLY' : 'g',
    'HIS' : 'h',
    'ILE' : 'i',
    'LEU' : 'l',
    'LYS' : 'k',
    'MET' : 'm',
    'PHE' : 'f',
    'PRO' : 'p',
    'SER' : 's',
    'THR' : 't',
    'TRP' : 'w',
    'TYR' : 'y',
    'VAL' : 'v',
}

"""
Task that takes a list of pdbs and processes the files extracting
geometry data from the files
"""
class ProcessPDBTask(Task):

    """
        Work function - expects a list of pdb file prefixes.
    """
    def _work(self, args):

        pdbs = args['pdbs']

        for data in pdbs:
            self.process_pdb(data)


    @transaction.commit_manually
    def process_pdb(self, data):
        try:

            code = data['code']
            filename = 'pdb%s.ent.Z' % code

            # update datastructure
            data['chains'] = []
            data['residues'] = {}

            # 1) parse with bioPython
            data = parseWithBioPython(filename, data)
            #print 'props: %s' % data

            # 2) parse with MMLib merging the dictionaries
            data = parseWithMMLib('%s/%s' % ('pdb', filename), data)

            # 3) Create/Get Protein and save values
            try:
                protein = ProteinModel.objects.get(code=code)
                print '  Existing protein: ', code
            except ProteinModel.DoesNotExist:
                print '  Creating protein: ', code
                protein = ProteinModel()
            protein.code       = code
            protein.threshold  = data['threshold']
            protein.resolution = data['resolution']
            protein.rfactor    = data['rfactor']
            protein.save()

            # 4) Get/Create Chains and save values
            chains = {}
            for chaincode in data['chains']:
                chainId = '%s%s' % (protein.code, chaincode)
                try:
                    chain = protein.chains.get(id=chainId)
                    print '   Existing Chain: %s' % chaincode
                except ChainModel.DoesNotExist:
                    print '   Creating Chain: %s' % chaincode
                    chain = ChainModel()
                    chain.id      = chainId
                    chain.protein = protein
                    chain.code    = chaincode
                    chain.save()

                    protein.chains.add(chain)
                #create dictionary of chains for quick access
                chains[chaincode] = chain


            # 5) iterate through residue data creating residues
            for id, residue_props in data['residues'].items():
                chain = chains[residue_props['chain_id']]

                # 5a) find the residue object so it can be updated or create a new one
                try:
                    residue = chain.residues.get(chainIndex=id)
                except ResidueModel.DoesNotExist:
                    #not found, create new residue
                    #print 'New Residue'
                    residue = ResidueModel()
                    residue.protein = protein
                    residue.chain   = chain
                    residue.chainID = chain.id[4]
                    residue.chainIndex = id
                    residue.oldID = id

                # 5b) copy properties into a residue object
                #     property keys should match property name in object
                for key, value in residue_props.items():
                    residue.__dict__[key] = value

                # 5c) save
                residue.save()
                chain.residues.add(residue)
        except Exception, e:
            import traceback
            exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
            print "*** print_tb:"
            traceback.print_tb(exceptionTraceback, limit=1, file=sys.stdout)

            print 'exception residue',e


        # 6) entire protein has been processed, commit transaction
        transaction.commit()


"""
Uncompress using the UNIX uncompress command.  The PDB files are stored 
used LZ algorithms i've tried will not decompress the files properly.
For now use the linux decompress
"""
def uncompress(file, src_dir, dest_dir):
    tempfile = '%s/%s' % (dest_dir, file)
    dest = '%s/%s' % (dest_dir, file[:-2])

    try:
        # copy the file to the tmp directory
        shutil.copyfile('%s/%s' % (src_dir,file), tempfile)

        # decompress using unix decompress.  
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


"""
Parse values from file that can be parsed using BioPython library
@return a dict containing the properties that were processed
"""
def parseWithBioPython(file, props):
    residues = props['residues']

    decompressedFile = None
    tmp = './tmp'
    pdb = './pdb'

    try:
        #create tmp workspace
        if os.path.exists(tmp):
            ownTempDir = False
        else:
            ownTempDir = True
            os.mkdir(tmp)

        #prep and open file
        decompressedFile = uncompress(file, pdb, tmp)

        if decompressedFile:

            structure = Bio.PDB.PDBParser().get_structure('pdbname', decompressedFile)

            # dssp can't do multiple models. if we ever need to, we'll have to 
            # iterate through them
            dssp = Bio.PDB.DSSP(model=structure[0], pdb_file=decompressedFile, dssp='dsspcmbi')

            #iterate residues
            res_old_id = None
            oldN        = None
            oldCA       = None
            oldC        = None

            for res in structure[0].get_residues():
                hetflag, res_id, icode = res.get_id()
                chain = res.get_parent().get_id()
                if hetflag != ' ':
                    oldN       = None
                    oldCA      = None
                    oldC       = None
                    continue

                # get properties using dssp
                residue, secondary_structure, accessibility, relative_accessibility = dssp[(chain, res_id)]

                # save in dictionary
                try:
                    residues[res_id]['ss'] = secondary_structure
                except:
                    # residue didn't exist yet
                    residues[res_id] = {'ss':secondary_structure}

                newN    = res['N'].get_vector()
                newCA   = res['CA'].get_vector()
                newC    = res['C'].get_vector()
                newCB   = res['CB'].get_vector() if res.has_id('CB') else None
                newO    = res['O'].get_vector()

                residues[res_id]['a1'] = NO_VALUE
                residues[res_id]['a2'] = NO_VALUE
                residues[res_id]['a3'] = NO_VALUE
                residues[res_id]['a4'] = NO_VALUE
                residues[res_id]['a5'] = NO_VALUE
                residues[res_id]['a6'] = NO_VALUE
                residues[res_id]['a7'] = NO_VALUE
                #residues[res_id]['L1'] = NO_VALUE
                residues[res_id]['psi'] = NO_VALUE
                residues[res_id]['ome'] = NO_VALUE
                residues[res_id]['phi'] = NO_VALUE

                if residues.has_key(res_old_id) and res_old_id+1 == res_id:
                    residues[res_id]['a6'] = math.degrees(calc_angle(oldCA,oldC,newN))
                    residues[res_id]['a7'] = math.degrees(calc_angle(oldO,oldC,newN))
                    residues[res_old_id]['psi']= math.degrees(calc_dihedral(oldN,oldCA,oldC,newN))
                    residues[res_old_id]['ome']    = math.degrees(calc_dihedral(oldCA,oldC,newN,newCA)) # or should it be on the newResObj?
                    residues[res_id]['a1']     = math.degrees(calc_angle(oldC,newN,newCA))
                    #residues[res_id]['L1']     = calc_angle(oldC,newN)
                    residues[res_id]['phi']    = math.degrees(calc_dihedral(oldC,newN,newCA,newC))

                #residues[res_id]['L2'] = calc_angle(newN,newCA)
                #residues[res_id]['L4'] = calc_angle(newCA,newC)
                #residues[res_id]['L5'] = calc_angle(newC,newO)
                residues[res_id]['a3'] = math.degrees(calc_angle(newN,newCA,newC))
                residues[res_id]['a5'] = math.degrees(calc_angle(newCA,newC,newO))

                if newCB:
                    residues[res_id]['a2'] = math.degrees(calc_angle(newN,newCA,newCB))
                    residues[res_id]['a4'] = math.degrees(calc_angle(newCB,newCA,newC))
                    #residues[res_id]['L3'] = calc_angle(newCA,newCB)

                res_old_id = res_id
                oldN       = newN
                oldCA      = newCA
                oldC       = newC
                oldO       = newO

    finally:
        #clean up any files in tmp directory no matter what
        if decompressedFile and os.path.exists(decompressedFile):
            os.remove(decompressedFile)

        if ownTempDir and os.path.exists(tmp):
            os.removedirs(tmp)

    return props


"""
Parse values from file that can be parsed using MMLib library
@return a dict containing the properties that were processed
"""
def parseWithMMLib(file, props):
    residues = props['residues']

    struct = mmLib.FileIO.LoadStructure(file=file)
    lowerCaseIndicators = ['H','G','E','T']

    props['code']       = struct.structure_id
    #props['threshold']  = struct.
    #props['resolution'] = struct.
    #props['rfactor']    = struct.

    for r in struct.iter_amino_acids():

        #get or create residue properties
        #this assumes residue property was already created in parseWithBioPython, Speed Optimization
        residueProps = residues[int(r.fragment_id)]

        # get offsets for residues left and right of this one
        next_res = r.get_offset_residue(1)
        if not next_res or int(next_res.fragment_id) != (int(r.fragment_id)+1): next_res = None
        prev_res = r.get_offset_residue(-1)
        if not prev_res or int(prev_res.fragment_id) != (int(r.fragment_id)-1): prev_res = None

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


        # determine if SS should be lowercase
        ss = residues[int(r.fragment_id)]['ss']

        if prev_res:
            prev_ss = residues[int(r.fragment_id)-1]['ss']
        else:
            prev_ss = None

        if next_res:
            next_ss = residues[int(r.fragment_id)+1]['ss']
        else:
            prev_ss = None

        if (prev_ss and next_ss and prev_ss.upper() != ss or ss != next_ss.upper()) and ss in lowerCaseIndicators:
            ss = ss.lower()

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

        #add chain id
        if not r.chain_id in props['chains']:
            props['chains'].append(r.chain_id)

        #add all properties to residue dict
        residueProps['aa']      = AA3to1[r.res_name]
        residueProps['id']      = r.fragment_id
        #residueProps['phi']     = math.degrees(r.phi)
        #residueProps['psi']     = math.degrees(r.psi)
        #residueProps['ome']     = math.degrees(r.ome)
        residueProps['L1']      = r.l1
        residueProps['L2']      = r.l2
        residueProps['L3']      = r.l3
        residueProps['L4']      = r.l4
        residueProps['L5']      = r.l5
        residueProps['L6']      = r.l6
        residueProps['L7']      = r.l7
        #residueProps['a1']      = math.degrees(r.a1)
        #residueProps['a2']      = math.degrees(r.a2)
        #residueProps['a3']      = math.degrees(r.a3)
        #residueProps['a4']      = math.degrees(r.a4)
        #residueProps['a5']      = math.degrees(r.a5)
        #residueProps['a6']      = math.degrees(r.a6)
        #residueProps['a7']      = math.degrees(r.a7)
        residueProps['ss']      = ss
        residueProps['chain_id']= r.chain_id
        residueProps['h_bond_energy'] = 0.00
        residueProps['zeta']    = math.degrees(r.zeta)
        residueProps['bg']      = r.b_gamma
        residueProps['bm']      = r.b_main
        residueProps['bs']      = r.b_side
        residueProps['zero']    = 0.00
        residueProps['chi']      = math.degrees(r.x1)

        #residueProps['chi2']      = math.degrees(r.x2)
        #residueProps['chi3']      = math.degrees(r.x3)
        #residueProps['chi4']      = math.degrees(r.x4)

    return props

"""
Initialize the dictionary for geometry data
"""
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


"""
Run if file is executed from the command line
"""
if __name__ == '__main__':
    task = ProcessPDBTask('Command Line Processor')

    args = ['1qe5']
    args = {'pdbs':[{'code':'153l','threshold':1, 'resolution':2,'rfactor':3}]}

    task._work(args)

