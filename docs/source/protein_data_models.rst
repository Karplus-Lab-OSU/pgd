*******************
Protein Data Models
*******************

.. .. image:: protein_model.png

Protein data model is composed of three classes: **Protein, Chain,** and **Residue**. The models are designed to represent the protein in the most compact way. With proper SQL indexes this is also the most efficient method of storing the data for search queries when the database contains greater than 2 million residues.

The models are defined as Django models meaning the data exists as both Python classes and SQL Tables. This allows the data to be accessed through SQL, or preferably as Python objects through the Django Query API. Django models can be used from any program provided your *Django environment is configured*

For more information on Django models:

    * `Django Model Reference
      <https://docs.djangoproject.com/en/dev/topics/db/models/#topics-db-models>`_
    * `Django Query Reference
      <https://docs.djangoproject.com/en/dev/topics/db/queries/#topics-db-queries>`_

-------
Protein
-------

Represents a Protein.

    * **code** - 4 letter code for protein as imported from pdb.
    * **rfactor**
    * **rfree**
    * **resolution**
    * **threshold**
    * **residues** - related field collection, returns queryset of all residues related to this protein
    * **chains** - related field colleciton, returns queryset of all chain related to this protein

-----
Chain
-----

Chain id's are stored to allow importing of multiple chains from a single protein. The search interface does not currently include selection of chains.

    * **protein** (protein_id) Foreign Key relation to the protein this residue belongs to. May be retrieved as a Protein object or identifier
    * **residues** - related field collection, returns queryset of all residues related to this protein

-------
Residue
-------

Represents a Residue (Amino Acid) belonging to a protein

    * **id** - unique identifier for residue, internal to database
    * **oldID** - identifier as listed in PDB file. May include icodes appended to the end.
    * **chainIndex** - numerical index of residue in the chain. All chain breaks are represented as a single gap in numbers. IE. 1,2,4,5.
    * **protein** (protein_id) Foreign Key relation to the protein this residue belongs to. May be retrieved as a Protein object or identifier
    * **prev** (prev_id) - Foreign Key relation to the previous residue in the chain, this may be retrieved as a Residue object or identifier
    * **next** (next_id) - Foreign Key relation to the next residue in the chain, this may be retrieved as a Residue object or identifier

    * **Bond Lengths**
        * L1 -
        * L2 -
        * L3 -
        * L4 -
        * L5 -

    * **Bond Angles**
        * A1 -
        * A2 -
        * A3 -
        * A4 -
        * A5 -
        * A6 -
        * A7 -

    * **Dihedral Angles**
        * Omega
        * Phi
        * Psi
        * Zeta

    * **Sidechain**
        * **X1 through X4** as defined by mmlib.

    * **B-Factor**
        * **Bm** - Average b-factor of mainchain atoms
        * **Bs** - Average b-factor of sidechain
        * **Bg** - Average b-factor of Cg atom if present
