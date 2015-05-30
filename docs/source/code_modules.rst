************
Code Modules
************

The Protein Geometry Database makes use of Django "apps". Apps are synonymous with a module or plugin. They are module bits of an application that are, for the most part, portable between other apps.

PGD is divided into three apps to allow code to be portable.

^^^^^^^^
PGD Core
^^^^^^^^

PGD Core defines the core data structures. There is no functionality contained within this app. It is intended to contain only the database so that it can be reused in other applications

^^^^^^^^^^
PGD Search
^^^^^^^^^^

PGD Search contains logic for searching and displaying data from PGD Core. This includes search logic, models, and the web front end.

^^^^^^^^^^^
PGD Splicer
^^^^^^^^^^^

PGD Splicer contains all of the code required for importing data from PDB files into objects defined within PGD Core
