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
    from pydra.config import load_settings
    pydra_settings = load_settings()



from datetime import datetime
import math
from math import sqrt
import os
import shutil
import sys
import time

import Bio.PDB
from Bio.PDB import calc_angle as pdb_calc_angle
from Bio.PDB import calc_dihedral as pdb_calc_dihedral
from django.db import transaction
from pydra.cluster.tasks.tasks import Task

from pgd.pgd_core.models import Protein as ProteinModel
from pgd.pgd_core.models import Chain as ChainModel
from pgd.pgd_core.models import Residue as ResidueModel
from pgd.pgd_core.models import Sidechain_ARG
from pgd.pgd_core.models import Sidechain_ASN
from pgd.pgd_core.models import Sidechain_ASP
from pgd.pgd_core.models import Sidechain_CYS
from pgd.pgd_core.models import Sidechain_GLN
from pgd.pgd_core.models import Sidechain_GLU
from pgd.pgd_core.models import Sidechain_HIS
from pgd.pgd_core.models import Sidechain_ILE
from pgd.pgd_core.models import Sidechain_LEU
from pgd.pgd_core.models import Sidechain_LYS
from pgd.pgd_core.models import Sidechain_MET
from pgd.pgd_core.models import Sidechain_PHE
from pgd.pgd_core.models import Sidechain_PRO
from pgd.pgd_core.models import Sidechain_SER
from pgd.pgd_core.models import Sidechain_THR
from pgd.pgd_core.models import Sidechain_TRP
from pgd.pgd_core.models import Sidechain_TYR
from pgd.pgd_core.models import Sidechain_VAL

from pgd.pgd_splicer.models import *
from pgd.pgd_splicer.chi import CHI_MAP
from pgd.pgd_splicer.sidechain import *


def NO_VALUE(field):
    """
    Helper function for determining the value to use when the field is an
    invalid value.
    """
    if field in ('bm','bg','bs'):
        return 0
    else:
        return None

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


aa_class = {
    'r':(Sidechain_ARG,'sidechain_ARG'),
    'n':(Sidechain_ASN,'sidechain_ASN'),
    'd':(Sidechain_ASP,'sidechain_ASP'),
    'c':(Sidechain_CYS,'sidechain_CYS'),
    'q':(Sidechain_GLN,'sidechain_GLN'),
    'e':(Sidechain_GLU,'sidechain_GLU'),
    'h':(Sidechain_HIS,'sidechain_HIS'),
    'i':(Sidechain_ILE,'sidechain_ILE'),
    'l':(Sidechain_LEU,'sidechain_LEU'),
    'k':(Sidechain_LYS,'sidechain_LYS'),
    'm':(Sidechain_MET,'sidechain_MET'),
    'f':(Sidechain_PHE,'sidechain_PHE'),
    'p':(Sidechain_PRO,'sidechain_PRO'),
    's':(Sidechain_SER,'sidechain_SER'),
    't':(Sidechain_THR,'sidechain_THR'),
    'w':(Sidechain_TRP,'sidechain_TRP'),
    'y':(Sidechain_TYR,'sidechain_TYR'),
    'v':(Sidechain_VAL,'sidechain_VAL')
}


class InvalidResidueException(Exception):
    """
    Exception identifying something wrong while processing
    a protein residue.  this is used to jump to a common error
    handling routine.
    """
    pass


class ProcessPDBTask(Task):
    """
    Task that takes a list of pdbs and processes the files extracting
    geometry data from the files.  The data is stored in Protein, Chain and
    Residue models and commited to the database.
    """
    total_proteins = 0
    finished_proteins = 0

    def progress(self):
        if not self.total_proteins:
            return 0
        return int((self.finished_proteins/float(self.total_proteins))*100)
            

    def work(self, **kwargs):
        """
        Work function - expects a list of pdb file prefixes.
        """
        # process a single protein dict, or a list of proteins
        pdbs = kwargs['data']
        if not isinstance(pdbs, list):
            pdbs = [pdbs]
        #print 'PDBS TO PROCESS:', pdbs
        self.total_proteins = len(pdbs)

        for data in pdbs:
            # only update pdbs if they are newer
            if self.pdb_file_is_newer(data):
                self.process_pdb(data)
            else:
                print 'INFO: Skipping up-to-date PDB: %s' % data['code']
            self.finished_proteins += 1

        print 'ProcessPDBTask - Processing Complete'

        # return only the code of proteins inserted or updated
        # we no longer need to pass any data as it is contained in the database
        # for now assume everything was updated
        codes = {'pdbs':[p['code'] for p in pdbs]}
        print codes
        return codes


    def pdb_file_is_newer(self, data):
        """
        Compares if the pdb file used as an input is newer than data already
        in the database.  This is used to prevent processing proteins
        if they do not need to be processed
        """
        code =  data['code']
        path = './pdb/pdb%s.ent.gz' % code.lower()
        print path
        if os.path.exists(path):
            pdb_date = datetime.fromtimestamp(os.path.getmtime(path))
            
        else:
            print 'ERROR - File not found'
            return False
        try:
            protein = ProteinModel.objects.get(code=code)
        except ProteinModel.DoesNotExist:
            # Protein not in database, pdb is new
            data['pdb_date'] = pdb_date
            return True

        data['pdb_date'] = pdb_date
        return protein.pdb_date < pdb_date


    @transaction.commit_manually
    def process_pdb(self, data):
        """
        Process an individual pdb file
        """
        try:
            residue_props = None
            code = data['code']
            chains_filter = data['chains'] if data.has_key('chains') else None
            filename = 'pdb%s.ent.gz' % code.lower()
            print '    Processing: ', code

            # update datastructure
            data['chains'] = {}

            # 1) parse with bioPython
            data = parseWithBioPython(filename, data, chains_filter)

            # 2) Create/Get Protein and save values
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
            protein.rfree      = float(data['rfree'])
            protein.pdb_date   = data['pdb_date']
            protein.save()

            # 3) Get/Create Chains and save values
            chains = {}
            for chaincode, residues in data['chains'].items():
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


                # 4) iterate through residue data creating residues
                comparison = lambda x,y: cmp(x['chainIndex'], y['chainIndex'])
                for residue_props in sorted(residues.values(), comparison):

                    # 4a) find the residue object so it can be updated or create a new one
                    try:
                        residue = chain.residues.get(oldID=str(residue_props['oldID']))
                    except ResidueModel.DoesNotExist:
                        #not found, create new residue
                        #print 'New Residue'
                        residue = ResidueModel()
                        residue.protein = protein
                        residue.chain   = chain
                        residue.chainID = chain.id[4]

                    # 4b) copy properties into a residue object
                    #     property keys should match property name in object
                    residue.__dict__.update(residue_props)

                    # 4c) set previous
                    if residue_props.has_key('prev'):
                        residue.prev = old_residue

                    # 4d) save
                    residue.save()
                    chain.residues.add(residue)
                    
                    # 4e) Update old_residue.next
                    if residue_props.has_key('prev'):
                        old_residue.next = residue
                        old_residue.save()

                    # 4f) find and update sidechain if needed
                    if 'sidechain' in residue_props:
                        klass, name = aa_class[residue_props['aa']]
                        try:
                            sidechain = getattr(residue, name)
                        except:
                            sidechain = klass()
                            sidechain.residue_id = residue.id
                        sidechain.__dict__.update(residue_props['sidechain'])
                        sidechain.save()
                    old_residue = residue
                print '    %s proteins' % len(residues)


        except Exception, e:
            import traceback
            exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
            print "*** print_tb:"
            print residue_props
            traceback.print_tb(exceptionTraceback, limit=10, file=sys.stdout)
            print 'EXCEPTION in Residue', code, e.__class__, e
            #self.logger.error('EXCEPTION in Residue: %s %s %s' % (code, e.__class__, e))
            print 'EXCEPTION in Residue: %s %s %s' % (code, e.__class__, e)
            #transaction.rollback()
            return

        # 5) entire protein has been processed, commit transaction
        transaction.commit()



def uncompress(file, src_dir, dest_dir):
    """
    Uncompress using the UNIX uncompress command.  PDB files are stored with
    GZIP.  Python supports GZIP but does not detect errors with incomplete
    files.  Its faster to just use the GZIP executable as it does not require
    loading the file into python and then dumping it back to a file.

    @param file: filename to decompress, does not include path
    @param src_dir: directory of file
    @param dest_dir: directory to write uncompressed file into
    """

    tempfile = '%s/%s' % (dest_dir, file)
    dest = '%s/%s' % (dest_dir, file[:-3])
    try:
        # copy the file to the tmp directory
        shutil.copyfile('%s/%s' % (src_dir,file), tempfile)

        # decompress using unix decompress.
        os.system('uncompress %s' % tempfile)

        # errors with uncompress won't be detected so we must
        # check for existence of the file
        if not os.path.exists(dest):
            raise Exception('File was not uncompressed')

    except Exception, e:
        print 'Exception while uncompressing file: %s - %s' % (tempfile, e)
        #clean up resulting file on errors
        if os.path.exists(dest):
            os.remove(dest)

        return False

    finally:
        # clean up temp file no matter what
        if os.path.exists(tempfile):
            print 'tempfile', tempfile
            os.remove(tempfile)

    return dest


def parseWithBioPython(file, props, chains_filter=None):
    """
    Parse values from file that can be parsed using BioPython library
    @return a dict containing the properties that were processed
    """
    chains = props['chains']

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

        if not decompressedFile:
            print 'ERROR: file not decompressed'
        else:

            structure = Bio.PDB.PDBParser().get_structure('pdbname', decompressedFile)

            # dssp can't do multiple models. if we ever need to, we'll have to 
            # iterate through them
            dssp = Bio.PDB.DSSP(model=structure[0], pdb_file=decompressedFile, dssp='dsspcmbi')

            for chain in structure[0]:
                chain_id = chain.get_id()

                # only process selected chains
                if chains_filter and not chain_id in chains_filter:
                    print 'Skipping Chain: %s' % chain_id
                    continue

                # construct structure for saving chain
                if not chain_id in props['chains']:
                    residues = {}
                    props['chains'][chain_id] = residues
                    print 'PROCESSING CHAIN [%s]' % chain, len(chain)

                newID = 0

                #iterate residues
                res_old_id = None
                oldN       = None
                oldCA      = None
                oldC       = None

                for res in chain:
                    try:
                        newID += 1
                        terminal = False
                        hetflag, res_id, icode = res.get_id()
                        
                        """
                        XXX Get the dictionary of atoms in the Main conformation.
                        BioPython should do this automatically, but it does not
                        always choose the main conformation.  Leading to some
                        Interesting results
                        """
                        atoms = {}
                        for atom in res.get_unpacked_list():
                            if atom.get_altloc() in ('A', ' '):
                                atoms[atom.name] = atom
                        
                        """
                        Exclude water residues
                        Exclude any Residues that are missing _ANY_ of the
                            mainchain atoms.  Any atom could be missing
                        """
                        all_mainchain = ('N' in atoms) and ('CA' in atoms) and ('C' in atoms) and ('O' in atoms)
                        if hetflag != ' ' or not all_mainchain:
                            raise InvalidResidueException('HetCode or Missing Atom')

                        """
                        Create dictionary structure and initialize all values.  All
                        Values are required.  Values that are not filled in will retain
                        the NO_VALUE value.

                        Store residue properties using OLD_ID as the key to ensure it is
                        unique.  We're including residues from all chains in the same
                        dictionary and chainindex may have duplicates.
                        """
                        old_id = res_id if icode == ' ' else '%s%s' % (res_id, icode)
                        try:
                            res_dict = residues[old_id]
                        except KeyError:
                            # residue didn't exist yet
                            res_dict = {}
                            residues[old_id] = res_dict
                            res_dict['oldID'] = old_id

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
                        try:
                            residue_dssp, secondary_structure, accessibility, relative_accessibility = dssp[(chain, (hetflag, res_id, icode)) ]
                        except KeyError, e:
                            import sys, traceback
                            t, v, tb = sys.exc_info()
                            traceback.print_tb(tb, limit=10, file=sys.stdout)
                            raise InvalidResidueException('KeyError in DSSP')

                        res_dict['chain_id'] = chain
                        res_dict['ss'] = secondary_structure
                        res_dict['aa'] = AA3to1[res.resname]
                        res_dict['h_bond_energy'] = 0.00

                        """
                        Get Vectors for mainchain atoms and calculate geometric angles,
                        dihedral angles, and lengths between them.
                        """
                        N    = atoms['N'].get_vector()
                        CA   = atoms['CA'].get_vector()
                        C    = atoms['C'].get_vector()
                        CB   = atoms['CB'].get_vector() if atoms.has_key('CB') else None
                        O    = atoms['O'].get_vector()

                        

                        if oldC:
                            # determine if there are missing residues by calculating
                            # the distance between the current residue and previous
                            # residue.  If the L1 distance is greater than 2.5 it 
                            # cannot possibly be the correct order of residues.
                            L1 = calc_distance(oldC,N)

                            if L1 < 2.5:
                                # properties that span residues
                                residues[res_old_id]['a6'] = calc_angle(oldCA,oldC,N)
                                residues[res_old_id]['a7'] = calc_angle(oldO,oldC,N)
                                residues[res_old_id]['psi'] = calc_dihedral(oldN,oldCA,oldC,N)
                                residues[res_old_id]['ome'] = calc_dihedral(oldCA,oldC,N,CA)
                                residues[res_old_id]['next'] = newID
                                res_dict['prev'] = residues[res_old_id]['chainIndex']
                                res_dict['a1']     = calc_angle(oldC,N,CA)
                                res_dict['phi']    = calc_dihedral(oldC,N,CA,C)
                                res_dict['L1'] = L1
                                terminal = False

                        if terminal:
                            # break in the chain, 
                            # 1) add terminal flags to both ends of the break so
                            #    the break can quickly be found.
                            # 2) skip a number in the new style index.  This allows
                            #    the break to be visible without checking the 
                            #    terminal flag
                            newID += 1
                            res_dict['terminal_flag'] = True
                            residues[res_old_id]['terminal_flag'] = True

                        # newID cannot be set until after we determine if it is terminal
                        res_dict['chainIndex'] = newID

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
                        for name in atoms:
                            if name in ('N', 'CA', 'C', 'O','OXT'):
                                main_chain.append(atoms[name].get_bfactor())
                            elif name in ('H'):
                                continue
                            else:
                                side_chain.append(atoms[name].get_bfactor())

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
                        sidechain = {}
                        calc_sidechain_lengths(res, sidechain)
                        calc_sidechain_angles(res, sidechain)
                        if sidechain:
                            res_dict['sidechain'] = sidechain
                        """
                        Reset for next pass.  We save some relationships which span two atoms.
                        """
                        res_old_id = old_id
                        oldN       = N
                        oldCA      = CA
                        oldC       = C
                        oldO       = O

                    except InvalidResidueException, e:
                        # something has gone wrong in the current residue
                        # indicating that it should be excluded from processing
                        # log a warning
                        print 'WARNING: Invalid residue - protein:%s  chain:%s   residue: %s  exception: %s' % (file, chain_id, res.get_id(), e)
                        if oldC:
                            residues[res_old_id]['terminal_flag'] = True
                            newID += 1
                        oldN       = None
                        oldCA      = None
                        oldC       = None
                        if residues.has_key(res_id):
                            del residues[res_id]

                print 'Processed %s residues' % len(residues)

    finally:
        #clean up any files in tmp directory no matter what
        if decompressedFile and os.path.exists(decompressedFile):
            os.remove(decompressedFile)

        if ownTempDir and os.path.exists(tmp):
            os.removedirs(tmp)

    return props


def initialize_geometry(residue, geometry_list, type):
    """
    Initialize the dictionary for geometry data
    """
    for item in geometry_list:
        if not residue.has_key(item) or residue[item] is None:
            if type == 'angle':
                residue[item] = NO_VALUE(item)
            elif type == 'length':
                residue[item] = NO_VALUE(item)
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


def calc_sidechain_lengths(residue, residue_dict):
    """
    Calculates Values for sidechain bond lengths. Uses a predefined list
    from sidechain.py, specifically bond_lengths.
    """
    try:
        mapping = bond_lengths[residue.resname]
        for i in range(len(mapping)):
            atom_names = mapping[i]
            try:
                sidechain_atoms = [residue[n].get_vector() for n in atom_names]
                sidechain_length = calc_distance(*sidechain_atoms)
                residue_dict['%s_%s' % (atom_names[0], atom_names[1])] = sidechain_length
            except KeyError:
                #missing an atom
                continue

    except KeyError:
        # this residue type does not have sidechain lengths
        pass


def calc_sidechain_angles(residue, residue_dict):
    """
    Calculates Values for sidechain bond angles. Uses a predefined list
    from sidechain.py, specifically bond_angles.
    """
    try:
        mapping = bond_angles[residue.resname]
        for i in range(len(mapping)):
            atom_names = mapping[i]
            try:
                sidechain_atoms = [residue[n].get_vector() for n in atom_names]
                sidechain_angle = calc_angle(*sidechain_atoms)
                residue_dict['%s_%s_%s' % (atom_names[0], atom_names[1], atom_names[2])] = sidechain_angle
            except KeyError:
                #missing an atom
                continue

    except KeyError:
        # this residue type does not have sidechain angles
        pass


if __name__ == '__main__':
    """
    Run if file is executed from the command line
    """
    #from pydra.cluster.worker import WorkerProxy
    import sys
    import logging
    import fileinput

    def process_args(args):
        return {'code':args[0],
                'chains':[c for c in args[1]],
                'threshold':float(args[2]),
                'resolution':float(args[3]),
                'rfactor':float(args[4]),
                'rfree':float(args[5])
                }

    task = ProcessPDBTask()
    
    logging.basicConfig(filename='ProcessPDB.log',level=logging.DEBUG)
    task.logger = logging
    #task.parent = WorkerProxy()

    pdbs = []

    argv = sys.argv
    if len(argv) == 1:
        print 'Usage:'
        print '   ProcessPDBTask code chains threshold resolution rfactor rfree [repeat]'
        print '       chains are a string of chain ids: ABCXYZ' 
        print ''
        print '   <cmd> | ProcessPDBTask --pipein'
        print '   piped protein values must be separated by newlines'
        sys.exit(0)
        
    elif len(argv) == 2 and argv[1] == '--pipein':
        for line in sys.stdin:
            pdbs.append(process_args(line.split(' ')))
            
    else:
        for i in range(1,len(argv),6):
            try:
                print argv[i:i+6]
                pdbs.append(process_args(argv[i:i+6]))
            except IndexError, e:
                print e
                print 'Usage: ProcessPDBTask.py code chain threshold resolution rfactor rfree...'
                sys.exit(0)
    
    task.work(**{'data':pdbs})