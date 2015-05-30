***************
Running Splicer
***************

Splicer is intended to be run on a `Pydra <https://code.osuosl.org/projects/pydra>`_ cluster. These are instructions and notes on running it and dealing with Pydra's immaturity.

Splicer can be `deployed <https://code.osuosl.org/projects/pydra>`_ and `run <https://code.osuosl.org/projects/pydra>`_ like any other task on a Pydra cluster.

Components of splicer can also be `run manually from the command line
<https://code.osuosl.org/projects/pgd/wiki/Designsplicercli>`_

---------------
Slow FTP Issues
---------------

The FTP server for PDB files is a very slow, rate limited, server located in the UK. PDB files are currently 1.8 gigabytes total for 16,000 proteins in PGD. It can take a long time to download this much data from the FTP server. This is handled in two ways:

-----------------------------------------
Maintaining Connections Between Workunits
-----------------------------------------

Each workunit is composed of downloading and processing a PDB file. Rather than disconnecting from the FTP server, connections are maintained until the last work unit is completed. This removes the overhead for connecting and disconnecting from the server

--------------------------
Only Downloading New Files
--------------------------

Checking dates is very fast, the MODTIME command completes almost instantly. This prevents uneeded downloading

--------------------------------
Storing files in a network share
--------------------------------

Pydra can't guarantee that future runs of Splicer will process the same set of proteins on the same hardware. This means that an up-to-date PDB could be mistaken for a PDB that doesn't exist. Storing the files on a shared filesystem ensures that regardless of which `Node <https://code.osuosl.org/projects/pydra#Node>`_ is assigned the workunit, it will find the same set of PDB files.

Note that this only matters when using

--------------------------
Workunit Thrashing Problem
--------------------------

There is an outstanding `bug in pydra that causes the node to crash when workunits complete too quickly <https://code.osuosl.org/projects/pydra#Node>`_. Splicer includes an option to batch process proteins to ensure that this does not happen. Eventually `batching workunits <https://code.osuosl.org/projects/pydra#Node>`_ will an automatic feature of Pydra

When running repeat runs of Pydra it is important to increase the workunit size to at least 500-1000. Because the date checks are very fast it will cycle through the existing proteins very quickly.

-----------------
Debugging Splicer
-----------------

Pydra `logs <https://code.osuosl.org/projects/pydra#Node>`_ most things that happen within it. A full `task history <https://code.osuosl.org/projects/pydra#Node>`_ can be viewed by clicking the history icon found on the pydra `tasks page <https://code.osuosl.org/projects/pydra#Node>`_. Clicking on a task instance gives you more details about the task including which workunits were successful and what their arguments were.

Workunits are logged individually and located in **/var/logs/pydra/archive**. The logs are aggregated from the Nodes after it is done with the entire task
