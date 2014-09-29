******************************
Optimization: In Memory Tables
******************************

Loading PGD tables in memory, if indexed properly, will greatly outperform on-disk tables. It is a solution dependent on a capable enough server.

----------------------
Indexing Memory Tables
----------------------

Memory tables will only outperform a properly indexed on-disk table if it is also indexed. Above a certain threshold a full table scan of a memory table is still slower than a binary search of an index, even when it is on disk.

Memory tables support indexing, but mysql does not correctly index when using btree. Btree is the default for on-disk tables, and is generally faster than a Hash index. however, mysql will generate an empty index when btree is used with memory tables. This means that any table joined with a memory table must also be in memory.

--------------------------
Parallelization of Queries
--------------------------

Memory addresses can be read by multiple threads simultaneously, unlike disks which require seeking back and forth. Without seek times to slow down simultaneous reads, multiple queries can be run on an in memory table at the same time. This can reduce the time required for a query to the longest in the set of queries

**Diagram**

Note that this becomes a limitation of CPU cores and the number of concurrent threads the server is capable of handling.

---------------------------------------
Startup and Django Configuration Issues
---------------------------------------

Memory tables do not persist through mysql restarts. They must be recreated and indexed every time the server starts. This needs to be automated so that when the server first starts it is able to check and create the tables if needed.

Django must also be told what this special table is. The choices are:

    * rename the table during the creation process. It might be impossible to determine the state of the tables though.
    * change the name of the table in the django configuration to match the in memory table. this can only be done prior to django loading, the django orm cannot be reinitialized

---------------
Growth Concerns
---------------

Growth of the database is a major concern when dealing with large memory tables. Ram is cheap, but not as cheap as disk space. The current PGD database requries about 1.5 gigabytes (1.1 for Residues) of ram to load the Protein and Residue table in memory. Two factors will increase growth:

    * Additional Fields - We have at least a 2 dozen new properties to add, which will add around 80% growth in the short term. There may be more fields added later also.
    * Additional Proteins added - Expected to be 10% growth per year.

The works out to the following projects:

    * current - 1.5 gigabytes
    * 6 Months - 2.4 gigabytes
    * 12 Months - 2.64 gigabytes
    * 2 Years - 2.9 gigabytes
    * 3 Years - 3.19 gigabytes
    * 4 Years - 3.51 gigabytes
    * 5 Years - 3.86 gigabytes

A server purchased now with 4 gigabytes of ram allocated for the just the memory table would last 5 years. This is two years past what is normally "end of life" for a server.

