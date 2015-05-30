*************************
Optimization: SQL Indexes
*************************

Like any database PGD relies on SQL Indexes for improved performance. This is a description of indexes used, and some that didn't work.

Please see the section on `Memory Table Indexes <https://code.osuosl.org/projects/pgd/wiki/Designmodelsoptimizationram#IndexingMemoryTables>`_ for more information about the types of indexes used with memory tables.

-------
Protein
-------

^^^^^^^^^^^^^^^^^
Primary Key Index
^^^^^^^^^^^^^^^^^

The primary key index is used when specific proteins have been selected by Code (primary key)

^^^^^^^^^^^^^^^^
Resolution Index
^^^^^^^^^^^^^^^^

The index on resolution is used in most cases. It filters large sets of proteins. The default query, with resolution <= 1.2 reduces the number of proteins from 16,000 to 2500.

As the number of proteins nears the total number of proteins MySQL will switch to performing a full table scan. Even with indexes on other fields it does not appear to use them.

^^^^^^^^^^^^^^
Failed Indexes
^^^^^^^^^^^^^^

We also attempted to create indexes with resolution and other fields. No noticeable increase was detected, MySQL always opted for the individual Resolution Index.

^^^^^^^^^^^^^^^^^^^^^^^^^
Protein Joined to Residue
^^^^^^^^^^^^^^^^^^^^^^^^^

When joining a Residue with its Protein an index on Residue.protein_id is used

^^^^^^^^^^^^^^
Failed Indexes
^^^^^^^^^^^^^^

We attenmpted to add additional fields to the protein_id index. It was actually slower than the protein_id index alone.

-------------------------
Residue Joined to Residue
-------------------------

Residues are joined to Residues for the previous and next relationships using the Primary Key index on Residue.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Note on Join Direction for Previous and Next
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Residues join from **residue_0.next** to **residue_1.id** ::

    SELECT * FROM pgd_core_residue r0 INNER JOIN pgd_core_residue r1 ON (r0.next = r1.id)

instead of **residue_0.id** to **residue_1.prev** ::

    SELECT * FROM pgd_core_residue r0 INNER JOIN pgd_core_residue r1 ON (r0.next = r1.id)

The latter appeared to be a faster query but is not possible with Django. The custom clause requires adding the where clause with `queryset.extra() <https://docs.djangoproject.com/en/dev/ref/models/querysets/#extra-select-none-where-none-params-none-tables-none-order-by-none-select-params-none>`_. But django will also add the original clause.
