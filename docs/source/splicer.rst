*******
Splicer
*******

Splicer is the tool used to import data from PDB files.

Splicer is built as a Parallel Pydra Task. Pydra is a framework for parallel and distributed jobs with python. Using Pydra allows splicer run times to be reduced nearly linearly across nodes in the cluster.

--------------
Task Structure
--------------

Splicer is made up of several subtasks. It is organized using both `ContainerTask <https://code.osuosl.org/projects/pydra#ContainerTask>`_ (sequential subtasks) and `ParallelTask <https://code.osuosl.org/projects/pydra#ParallelTask>`_ to divide work up and parallelize it
It is organized like so:

    * Splicer (!ContainerTask)

        * Selector
        * ProcessProteinTask(!ParallelTask)

            * FTPUpdate
            * Processor
            * SegmentBuilder

^^^^^^^^
Selector
^^^^^^^^

Selector downloads lists of pdbs from a source, currently `Dunbracks PICSES <http://dunbrack.fccc.edu/PISCES.php>`_, and processes the list into a python object.

There are two selectors:
    * DunbracksSelector - selects records from dunbracks culled lists
    * PDBSelectSelector - selects records from pdbselect lists. This was replaced with dunbrack's selector

Both selectors filter on and parse protein level properties from the list.

^^^^^^^^^^
FTP Update
^^^^^^^^^^

Synchronizes a local cache of PDBs with a remote ftp server. By default it will synchronize all files located on the remote directory, optionally it can be given a list of PDB codes to synchronize.

FTP Update compares dates to the milisecond. This requires using the MODTIME (mdtm) FTP command.

^^^^^^^^^
Processor
^^^^^^^^^

The processor parses a pdb file and extracts properties from it. The current implementation uses *BioPython* library and DSSP. This component is `covered in depth <https://code.osuosl.org/projects/pgd/wiki/Designsplicerprocessor>`_.

^^^^^^^^^^^^^^^
Segment Builder
^^^^^^^^^^^^^^^

This is a defunct subtask that was used to build the *de-normalized segment table*
