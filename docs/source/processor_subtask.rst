**********************
Splicer Processor Task
**********************

The **Splicer Processor Task**, also known as **ProcessPDBTask** accepts a PDB file and processes it into `PGD Protein Models <https://code.osuosl.org/projects/pgd/wiki/Designmodels>`_.

=============================
Running From the command line
=============================

Like all splicer components *ProcessPDBTask* can be run from the commandline to simplify debugging the task. It requires the follow arguments:
    * PDB Code - 4 character alphanumeric code for a protein.
    * Threshold - Float
    * Resolution - Float
    * Rfactor - Float
    * Rfree - Float

When debugging only the PDB code need be a real value. The other values are required, but not validated.

=========
Libraries
=========

    * `BioPython <http://biopython.org/wiki/Main_Page>`_ - A library that can parse PDB files and contains various functions within its `PDB API <http://www.biopython.org/DIST/docs/api/>`_ for extracting data from them.
    * DSSP - A program for calculating secondary structure. BioPython has bindings for this program.

============
Parsing PDBs
============

Some properties are available as simple properties using the `Residue <http://www.biopython.org/DIST/docs/api/Bio.PDB.Residue.Residue-class.html>`_ class within BioPython. Most others require calculations involving individual `Atoms <http://www.biopython.org/DIST/docs/api/Bio.PDB.Atom.Atom-class.html>`_ within a Residue

----------------------------
Parsing Geometric Properties
----------------------------

Geometric properties must be calculated from raw atom data. BioPython supplies several functions for calculating functions between atom vectors.
    * **calc_length(vector, vector)** - Calculates distance between atoms in 3D space (supplied by PGD)
    * `calc_angle <http://www.biopython.org/DIST/docs/api/Bio.PDB.Vector%27-module.html#calc_angle>`_ - Calculates the angle between 3 atoms.
    * `calc_dihedral <http://www.biopython.org/DIST/docs/api/Bio.PDB.Vector%27-module.html#calc_dihedral>`_ - Calculates the dihedral (torsion) angle between 4 atoms

These functions require vectors which can be retrieved using `Atom.get_vector() <http://www.biopython.org/DIST/docs/api/Bio.PDB.Atom.Atom-class.html#get_vector>`_

==========
Example a3
==========

::
    N  = residue['N'].get_vector()
    CA = residue['CA'].get_vector()
    C  = residue['C'].get_vector()

    a3 = calc_angle(N,CA,C)

============
Example: Ome
============

::
    oldCA = prev_residue['CA'].get_vector()
    oldC  = prev_residue['C'].get_vector()
    N  = residue['N'].get_vector()
    CA = residue['CA'].get_vector()

    ome = calc_dihedral(oldCA,oldC,N,CA)

---------------------------
Parsing Averaged Properties
---------------------------

Several properties are presented as averages across the **Main-Chain, *Side-Chain,** and **Carbon Gamma** atoms.

    * **Main-Chain** - atoms N, C-alpha, C, O, OXT
    * **Side-Chain** - all other atoms excluding Main-Chain, C-gamma, and HETs (water)
    * **Carbon-gamma** - also known as C-gamma or Cg. A single atom.

This properties are calculated as the average of an atom level property across the atoms in this group.

============================
Example B-factor: Bm, Bs, Bg
============================

::

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

-----------------------------
Parsing Side Chain Properties
-----------------------------

Side chains are different for each type of atom. They require a map of connections to determine which atoms require angles, lengths, and dihedral angles calculated.

Currently **Chi1, Chi2, Chi3, and Chi4*are calculated using a map in *chi.py**. Sidechain lengths and angles will be added later, requiring an additional map of connections between atoms.

===============
Update Checking
===============

PDB Processor checks for updates when processing proteins. The `Protein Model <https://code.osuosl.org/projects/pgd/wiki/Designmodels>`_ contains a timestamp which corresponds to the update timestamp on the PDB file it was imported from. A protein will only be processed if the PDB file is newer.
