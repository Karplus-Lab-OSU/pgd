**************
Importing Data
**************

-------------------------------------
Running Splicer From The Command Line
-------------------------------------

Splicer can be run from the command line. It requires that several steps be run separately.

All commands should be run from the project root (directory with settings.py in it).

^^^^^^^^^^^^^^^^^^
Selecting Proteins
^^^^^^^^^^^^^^^^^^

Proteins must first be selected. Default filtering settings will be used for threshold, resolution and r_factor. ::

    1 ./pgd_splicer/dunbrack_selector.py

This will return information about the the parameters used, the files proteins were selected from, and a list of proteins in the following format: ::

    code chains threshold resolution rfactor rfree

Save the selection into a file::

    1 ./pgd_splicer/dunbrack_selector.py --pipeout > selection.txt

^^^^^^-
Options
^^^^^^-

    * --pipeout - will only output the data. This should be used if you would like to create output suitable for input into one of the other steps

^^^^^^^^^^^^^^^^^^^^^
Downloading PDB Files
^^^^^^^^^^^^^^^^^^^^^

PDB files are downloaded from an FTP site using **ftpupdate.py**. This script will synchronize **./pdb** with the remote FTP server. Only new files will be downloaded, but it will check the timestamps on all files.

This is a time-consuming step. Be prepared to wait for approximately two days for this to complete on a fresh local copy, or one day on an update. ::

    1 ./pgd_splicer/ftpupdate.py code [code...]

To only grab the proteins which are selected (and cut down massively on consumed bandwidth and time), try::

    1./pgd_splicer/ftpupdate.py --pipein < selection.txt

^^^^^^^^^^^^^^^^^^^^
Processing PDB Files
^^^^^^^^^^^^^^^^^^^^

PDB files can be imported into the database with **ProcessPDBTask.py**. Multiple proteins can be fed as commands to be imported. Errors will be written to **ProcessPDB.log**. ::

    1 ./pgd_splicer/ProcessPDBTask.py code chains threshold resolution rfactor rfree [repeat]

As before, a selection can be piped in::

    1 ./pgd_splicer/ProcessPDBTask.py --pipein < selection.txt

Expect this to take a few days as well.

^^^^^^^^^^
Parameters
^^^^^^^^^^

Parameters are all required, and may be repeated for multiple proteins.

    * **code** - protein code to import, should be all uppercase
    * **chains** - list of chains to import, should be a string of chain ids. (ie. ABCDEF). The string should not have quotes around it.
    * **threshold, resolution, rfactor, rfree** - the value for these fields. These properties are retrieved from the selection script so they are included as input for processing the protein.

^^^^^^^
Options
^^^^^^^

    * --pipein - input will be read from a pipe instead of arguments. proteins in the list should be separated by newlines.

-------
Example
-------

Some examples. Intermediate output is saved to a text file so that it can be examined later.

^^^^^^^^^^^
Full Import
^^^^^^^^^^^

Update all proteins regardless of whether the file was downloaded by **ftpupdate**. *ProcessPDBTask* will still check the update date and exclude pdbs that are not new. ::

    1 ./pgd_splicer/dunbrack_selector.py --pipeout > selected_proteins.txt
    2 ./pgd_splicer/ftpupdate.py --pipein < selected_proteins.txt
    3 ./pgd_splicer/ProcessPDBTask.py --pipein < selected_proteins.txt

^^^^^^^^^^^^^^^
Update Only New
^^^^^^^^^^^^^^^

Update only proteins for which we have a new FTP file. ::

    1 ./pgd_splicer/dunbrack_selector.py --pipeout > selected_proteins.txt
    2 ./pgd_splicer/ftpupdate.py --pipein --pipeout < selected_proteins.txt > updated_proteins.txt
    3 ./pgd_splicer/ProcessPDBTask.py --pipein < updated_proteins.txt

^^^^^^^^^^^^^^^^^^^^^^^^
Update Skipping Download
^^^^^^^^^^^^^^^^^^^^^^^^

If all the pdb files are already downloaded you may skip the FTP step to save time. ::

    1 ./pgd_splicer/dunbrack_selector.py --pipeout > selected_proteins.txt
    2 ./pgd_splicer/ProcessPDBTask.py --pipein < selected_proteins.txt

or as a single command::

    1 ./pgd_splicer/dunbrack_selector.py --pipeout | ./pgd_splicer/ProcessPDBTask.py --pipein

