******
Browse
******

Browse displays all residues for each segment in a result set. Results a

--------------
Selecting Data
--------------

Queries using Django's ORM focus on a single object. Accessing related fields such as **Residue.prev** or **Residue.next** result in a second query to resolve those objects. This means to display a single segment of length 5 you must do 4 additional queries. To display 25 segments per page would require 125 additional queries.

The solution is to calculate the list of residues for each segment and select
them by range. For example, to retrieve a segment of length 3 where i-3 select all residues from i-1 to i+1 (2 through 4). ::

    Select * from pgd_core_residue where chainIndex between 2 and 4

This allows all residues in a segment to be retrieved in a single query. This
reduces the number of queries to the number of residues in the result set. This may still be a large number of queries but each page is limited to 25 segments per page.
