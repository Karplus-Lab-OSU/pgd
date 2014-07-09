*****************
Search Statistics
*****************

The statistics page provides statistics of properties across all results for a specific residue index. Results are grouped by Amino Acid Type, and or Secondary Structure Type.

-------
Queries
-------

The statistics page makes use of aggregate functions to calculate values. Because of the different types of grouping and statistics, it requires several queries to retrieve all of the statistics.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Secondary Structure Counts and Amino Acid Totals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This query produces three values:

    * Counts of residues per Secondary Structure Type, per Amino Acid Type.
    * Counts of residues per Amino Acid Type.
    * Total count of all residues in search.

The SQL required is ::

    Select count(*) as count, aa_type, ss_type from pgd_core_residue GROUP BY aa_type, ss_type WITH ROLLUP;

The **WITH ROLLUP** clause instructs mysql to include totals for the fields. This produces the total per amino acid type, and total count of all residues. WITH ROLLUP only affects the first field in the GROUP BY class so totals per Secondary Structure Type are not produced by this query.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^--
Secondary Structure Total Counts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^--

Count of residues per Secondary Structure Type. ::

    Select count(*) as count, ss_type from pgd_core_residue GROUP BY ss;

^^^^^^^^^^^^^^^^
Field Statistics
^^^^^^^^^^^^^^^^

Min, Max, Average, and Standard Deviation are calculated for every residue, per Amino Acid Type, include totals. ::

    SELECT MIN(a1), MAX(a1), AVG(a1), STDDEV(a1) from pgd_core_residue GROUP_BY aa WITH ROLLUP;

^^^^^^^^^^^^^^^^^^^^^^^^^^
Dihedral Angles Statistics
^^^^^^^^^^^^^^^^^^^^^^^^^^

Dihedral angles (ome, zeta, phi, psi, etc.) require the use of `Special Statistics for Dihedral Angles <https://code.osuosl.org/projects/pgd/wiki/Designmodelsoptimizationsql_aggregates#StatisticsforDihedralAngles>`_. Average is calculated in the field statistics query, but standard deviation requires a second query.

------------
Optimization
------------

Statistics uses the following optimizations:

    * `SQL Aggregate Functions <https://code.osuosl.org/projects/pgd/wiki/Designmodelsoptimizationsql_aggregates>`_ to reduce network transport overhead
    * `Parallel Queries <https://code.osuosl.org/projects/pgd/wiki/Designmodelsoptimizationram#ParallelizationofQueries>`_ to run calculations simultaneously
