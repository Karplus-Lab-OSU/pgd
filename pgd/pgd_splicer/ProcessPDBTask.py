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

from pydra_server.cluster.tasks import Task

from pgd_splicer.models import *
from pgd_core.models import Protein as ProteinModel
from pgd_core.models import Chain as ChainModel
from pgd_core.models import Residue as ResidueModel
from chi import CHI_MAP

from django.db import transaction

import math, sys, mmLib.Structure, mmLib.FileIO, mmLib.AtomMath
from math import sqrt
import Bio.PDB
from Bio.PDB import calc_angle as pdb_calc_angle
from Bio.PDB import calc_dihedral as pdb_calc_dihedral
import shutil
import os

import logging
logger = logging.getLogger('root')

class ICodeException(Exception):
    pass

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


class ProcessPDBTask(Task):
    """
    Task that takes a list of pdbs and processes the files extracting
    geometry data from the files.  The data is stored in Protein, Chain and
    Residue models and commited to the database.
    """

    def _work(self, **kwargs):
        """
        Work function - expects a list of pdb file prefixes.
        """

        # process a single protein dict, or a list of proteins
        pdbs = kwargs['data']
        if not isinstance(pdbs, list):
            pdbs = [pdbs]
        print 'PDBS TO PROCESS:', pdbs

        for data in pdbs:
            self.process_pdb(data)

        # return only the code of proteins inserted or updated
        # we no longer need to pass any data as it is contained in the database
        # for now assume everything was updated
        return {'pdbs':[p['code'] for p in pdbs]}

    @transaction.commit_manually
    def process_pdb(self, data):
        """
        Process an individual pdb file
        """
        try:
            code = data['code']
            print 'DATA', data
            filename = 'pdb%s.ent.gz' % code.lower()
            print '    Processing: ', code

            # update datastructure
            data['chains'] = []
            data['residues'] = {}

            # 1) parse with bioPython
            try:
                data = parseWithBioPython(filename, data)
            except ICodeException:
                logger.warn('PDB File had insertion codes: %s' % code)
                return
            #print 'props: %s' % data

            # 2) parse with MMLib merging the dictionaries
            #data = parseWithMMLib('%s/%s' % ('pdb', filename), data)

            # 3) Create/Get Protein and save values
            try:
                protein = ProteinModel.objects.get(code=code)
                print '  Existing protein: ', code
            except ProteinModel.DoesNotExist:
                print '  Creating protein: ', code
                protein = ProteinModel()
            protein.code       = code
            protein.threshold  = float(data['threshold'])
            protein.resolution = float(data['resolution'])
            protein.rfactor    = float(data['rfactor'])
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
            traceback.print_tb(exceptionTraceback, limit=10, file=sys.stdout)

            logger.error('EXCEPTION in Residue', code, e)
            transaction.rollback()

        # 6) entire protein has been processed, commit transaction
        transaction.commit()


"""
Uncompress using the UNIX uncompress command.  The PDB files are stored 
used LZ algorithms i've tried will not decompress the files properly.
For now use the linux decompress
"""
def uncompress(file, src_dir, dest_dir):
    tempfile = '%s/%s' % (dest_dir, file)
    dest = '%s/%s' % (dest_dir, file[:-3])
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


def parseWithBioPython(file, props):
    """
    Parse values from file that can be parsed using BioPython library
    @return a dict containing the properties that were processed
    """
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

            for chain in structure[0]:
                if not chain.get_id() in props['chains']:
                    props['chains'].append(chain.get_id())
                    print 'CHAIN [%s]' % chain

            for res in structure[0].get_residues():
                hetflag, res_id, icode = res.get_id()

                if icode != ' ':
                    raise ICodeException()

                if hetflag != ' ':
                    oldN       = None
                    oldCA      = None
                    oldC       = None
                    continue

                #print 'RES:', res_id, icode, res.parent

                """
                Create dictionary structure and initialize all values.  All
                Values are required.  Values that are not filled in will retain
                the NO_VALUE value.
                """
                try:
                    res_dict = residues[res_id]
                except KeyError:
                    # residue didn't exist yet
                    res_dict = {}
                    residues[res_id] = res_dict
                    res_dict['id'] = res_id

                length_list = ['L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7','bg','bs','bm']
                angles_list = ['a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7']
                dihedral_list = ['psi', 'ome', 'phi', 'zeta','chi1','chi2','chi3','chi4']
                initialize_geometry(res_dict, length_list, 'length')
                initialize_geometry(res_dict, angles_list, 'angle')
                initialize_geometry(res_dict, dihedral_list, 'angle')


                """
                Get Properties from DSSP and other per residue properties
                """
                chain = res.get_parent().get_id()
                residue_dssp, secondary_structure, accessibility, relative_accessibility = dssp[(chain, res_id)]
                res_dict['chain_id'] = chain
                res_dict['ss'] = secondary_structure
                res_dict['aa'] = AA3to1[res.resname]
                res_dict['h_bond_energy'] = 0.00

                """
                Get Vectors for mainchain atoms and calculate geometric angles,
                dihedral angles, and lengths between them.
                """
                N    = res['N'].get_vector()
                CA   = res['CA'].get_vector()
                C    = res['C'].get_vector()
                CB   = res['CB'].get_vector() if res.has_id('CB') else None
                O    = res['O'].get_vector()

                if residues.has_key(res_old_id) and res_old_id+1 == res_id:
                    # properties that span residues
                    res_dict['a6'] = calc_angle(oldCA,oldC,N)
                    res_dict['a7'] = calc_angle(oldO,oldC,N)
                    residues[res_old_id]['psi'] = calc_dihedral(oldN,oldCA,oldC,N)
                    residues[res_old_id]['ome'] = calc_dihedral(oldCA,oldC,N,CA)
                    res_dict['a1']     = calc_angle(oldC,N,CA)
                    res_dict['L1']     = calc_distance(oldC,N)
                    res_dict['phi']    = calc_dihedral(oldC,N,CA,C)

                res_dict['L2'] = calc_distance(N,CA)
                res_dict['L4'] = calc_distance(CA,C)
                res_dict['L5'] = calc_distance(C,O)
                res_dict['a3'] = calc_angle(N,CA,C)
                res_dict['a5'] = calc_angle(CA,C,O)

                if CB:
                    res_dict['a2'] = calc_angle(N,CA,CB)
                    res_dict['a4'] = calc_angle(CB,CA,C)
                    res_dict['L3'] = calc_distance(CA,CB)
                    res_dict['zeta'] = calc_dihedral(CA, N, C, CB)


                """
                Calculate Bg - bfactor of the 4th atom in Chi1.
                """
                try:
                    atom_name = CHI_MAP[res.resname][0][3]
                    res_dict['bg'] = res[atom_name].get_bfactor()
                except KeyError:
                    # not all residues have chi
                    pass


                """
                Other B Averages
                    Bm - Average of bfactors in main chain.
                    Bm - Average of bfactors in side chain.
                """
                main_chain = []
                side_chain = []
                for a in res.child_list:
                    if a.name in ('N', 'CA', 'C', 'O','OXT'):
                        main_chain.append(a.get_bfactor())
                    elif a.name in ('H'):
                        continue
                    else:
                        side_chain.append(a.get_bfactor())

                if main_chain != []:
                    res_dict['bm'] = sum(main_chain)/len(main_chain)

                if side_chain != []:
                    res_dict['bs'] = sum(side_chain)/len(side_chain)


                """
                Calculate CHI values.  The mappings for per peptide chi's are stored
                in a separate file and a function is used to calculate the chi based
                based on the peptide of this residue and the lists of atoms in the
                chi mappings.
                """
                calc_chi(res, res_dict)
                if res_dict['chi1']:
                    res_dict['chi'] = res_dict['chi1']


                """
                Reset for next pass.  We save some relationships which span two atoms.
                """
                res_old_id = res_id
                oldN       = N
                oldCA      = CA
                oldC       = C
                oldO       = O

    finally:
        #clean up any files in tmp directory no matter what
        if decompressedFile and os.path.exists(decompressedFile):
            os.remove(decompressedFile)

        if ownTempDir and os.path.exists(tmp):
            os.removedirs(tmp)

    return props


def parseWithMMLib(file, props):
    """
    Parse values from file that can be parsed using MMLib library
    @return a dict containing the properties that were processed
    """
    struct = mmLib.FileIO.LoadStructure(file=file)
    for r in struct.iter_amino_acids():

        #get or create residue properties
        #this assumes residue property was already created in parseWithBioPython, Speed Optimization
        try:
            residueProps = residues[int(r.fragment_id)]
        except Exception,e:
            print '>>>>>>>>>>>>>>>>', r.chain_id, file
            raise e

    return props


def initialize_geometry(residue, geometry_list, type):
    """
    Initialize the dictionary for geometry data
    """
    for item in geometry_list:
        if not residue.has_key(item) or residue[item] is None:
            if type == 'angle':
                residue[item] = math.degrees(math.radians(NO_VALUE))
            elif type == 'length':
                residue[item] = NO_VALUE
            else:
                print "Don't know how to deal with type", type


def calc_distance(atom1, atom2):
    """
    Calculates distance between atoms because this is not built into BIOPython

    scribed from http://www.scribd.com/doc/9816032/BioPython-for-Bioinfo
    """
    dx = atom1[0] - atom2[0]
    dy = atom1[1] - atom2[1]
    dz = atom1[2] - atom2[2]
    return sqrt(dx*dx + dy*dy + dz*dz)


def calc_angle(atom1, atom2, atom3):
    """
    overridding pdb version of function to return values converted to degress
    """
    return math.degrees(pdb_calc_angle(atom1, atom2, atom3))


def calc_dihedral(atom1, atom2, atom3, atom4):
    """
    overridding pdb version of function to return values converted to degress
    """
    return math.degrees(pdb_calc_dihedral(atom1, atom2, atom3, atom4))


def calc_chi(residue, residue_dict):
    """
    Calculates Values for CHI using the predefined list of CHI angles in
    the CHI_MAP.  CHI_MAP contains the list of all peptides and the atoms
    that make up their different chi values.  This function will process
    the values known to exist, it will also skip chi values if some of the
    atoms are missing
    """
    try:
        mapping = CHI_MAP[residue.resname]
        for i in range(len(mapping)):
            chi_atom_names= mapping[i]
            try:
                chi_atoms = [residue[n].get_vector() for n in chi_atom_names]
                chi = calc_dihedral(*chi_atoms)
                residue_dict['chi%i'%(i+1)] = chi
            except KeyError:
                #missing an atom
                continue

    except KeyError:
        # this residue type does not have chi
        pass


if __name__ == '__main__':
    """
    Run if file is executed from the command line
    """
    task = ProcessPDBTask('Command Line Processor')

    args = {'pdbs':[{'code':'153l','threshold':1, 'resolution':2,'rfactor':3}]}

    task._work(args)

