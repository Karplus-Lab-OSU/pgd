****************************************************
Attempted Optimization: De-normalizing Residue Table
****************************************************

Searching segments of residues requires joining the residue table on itself, numerous times. Even when indexed properly joins can be slow. To remove the need for joins we de-normalized **pgd_core_residue** into **pgd_search_segment**. This table contained all properties, for each residue, for each possible segment in a protein.

This optimization greatly out-performed joins with less than 20,000 results. However, when the result set increased to 300,000 or more results it was twice as slow. The reason was that the table was too large to perform an entire table scan. It only was sped up when there was an index touching every field a search had a clause for. It was not feasible to build an index touch all 200+ fields to allow quick searching.

----------
Table Size
----------

The table size was not completely unmanageable, it was only 7 gigabytes. This was fine for on-disk, but would not scale well in memory.
