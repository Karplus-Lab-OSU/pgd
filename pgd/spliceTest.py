#!/usr/bin/python2.5
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

from Bio.PDB import *
from pgd_core.models import *

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

AAchis =  {
    'ALA' : lambda x: (),
    'ARG' : 'r',
    'ASN' : 'n',
    'ASP' : 'd',
    'CYS' : 'c',
    'GLU' : 'e',
    'GLN' : 'q',
    'GLY' : lambda x: (),
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
    'VAL' : lambda x: (
                calc_dihedral(x[''].get_vector(),x[''],x[''],x['']),
            ),
}

filename = 'pdb/pdb1fxk.ent'
prot = PDBParser().get_structure('X',filename)
dssp = DSP(model=model, pdb_file=filename, dssp='dsspcmbi')

protObj = Protein(
    code        = prot.id,
    threshold   = 1, #?
    resolution  = 1.5, #?
    rfactor     = 1.5, #?
)
protObj.save()

model = prot[0]
for chain in model:

    chainObj = Chain(
        id      = protObj.id+chain.id,
        protein = protObj,
        code    = chain.id,
    )
    chainObj.save()

    oldResObj   = None
    oldResId    = None
    oldN        = None
    oldCA       = None
    oldC        = None

    for newRes in chain:

        resname = newRes.get_resname()

        if resname.has_key(resname):
            newResId    = newRes.id[1]
            newResObj   = Residue(
                protein     = protObj,
                chain       = chainObj,
                aa          = AA3to1[resName],
                chainID     = chainObj.id,
                chainIndex  = newResID,
                ss          = dssp[(chain.id, newResId)],
                bm          = sum()/4.0
                # now add the b* stuff!
            )

            newN    = newRes['N'].get_vector()
            newCA   = newRes['CA'].get_vector()
            newC    = newRes['C'].get_vector()
            newCB   = newRes['CB'].get_vector() if newRes.has_id('CB') else None
            newO    = newRes['O'].get_vector()

            if oldResId != None and OldResId+1 == newResId:
                oldResObj.a6    = calc_angle(oldCA,oldC,newN)
                oldResObj.a7    = calc_angle(oldO,oldC,newN)
                oldResObj.psi   = calc_dihedral(oldN,oldCA,oldC,newN)
                newResObj.ome   = calc_dihedral(oldCA,oldC,newN,newNCA) # or should it be on the newResObj?
                newResObj.a1    = calc_angle(oldC,newN,newCA)
                newResObj.L1    = calc_angle(oldC,newN)
                newResObj.phi   = calc_dihedral(oldC,newN,newCA,newC)

            if oldResObj:
                oldResObj.save()

            newResObj.L2 = calc_angle(newN,newCA)
            newResObj.L4 = calc_angle(newCA,newC)
            newResObj.L5 = calc_angle(newC,newO)
            newResObj.a3 = calc_angle(newN,newCA,newC)
            newResObj.a5 = calc_angle(newCA,newC,newO)

            if newCB:
                newResObj.a2 = calc_angle(newN,newCA,newCB)
                newResObj.a4 = calc_angle(newCB,newCA,newC)
                newResObj.L3 = calc_angle(newCA,newCB)

            if newRes.has_id('CG'):
                newResObj.bg = newRes['CG'].get_bfactor() #double check this?
            newResObj.bs  =   (
                newN.get_bfactor() +
                newCA.get_bfactor() +
                newC.get_bfactor() +
                newO.get_bfactor()
            )/4

            sideSet = [a.get_bfactor() for a in newRes.child_list if (
                a != newN and a != newCA and a != newC and a != newO
            )]

            if sideSet != []:
                bs      = sum(sideSet)/len(sideSet)

            oldResObj   = newResObj
            oldResId    = oldResId
            oldN        = newN
            oldCA       = newCA
            oldC        = oldC

        elif oldResId:
            oldResObj   = None
            oldResId    = None
            oldN        = None
            oldCA       = None
            oldC        = None
    else:
        if oldResObj:
            oldResObj.save()
