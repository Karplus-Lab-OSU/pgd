*********
Data Dump
*********

Data dump is a dump of all fields of all segments returned by a search.

--------------
Selecting Data
--------------

Queries using Django's ORM focus on a single object. Accessing related fields such as **Residue.prev** or **Residue.next** result in a second query to resolve those objects. This means to display a single segment of length 5 you must do 4 additional queries.

The solution is to calculate the list of residues for each segment and select them by range. For example, to retrieve a segment of length 3 where i-3 select all residues from i-1 to i+1 (2 through 4). ::

    Select * from pgd_core_residue where chainIndex between 2 and 4

This allows all residues in a segment to be retrieved in a single query. This reduces the number of queries to the number of residues in the result set. This may still be a large number of queries but they may be run in parallel and the buffered nature of datadump spreads them out.

-----------------
Buffered Response
-----------------

Datadump uses an iterable class which uses threads to buffer data. The buffering threads use a paginator to split a result set into pieces. This allows downloading to start almost immediately after the button is clicked, rather than waiting for the entire dump to be written to memory.

This is not a perfect solution. Python uses green threads so threads are not truly concurrent. Threads will often be starved and may simply alternate between filling the buffer and emptying it. This is mainly an issue with the threads fighting over the lock. There may be a better solution using double-buffering.
