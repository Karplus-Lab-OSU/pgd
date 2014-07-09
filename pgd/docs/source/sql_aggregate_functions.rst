**************
SQL Aggregates
**************

SQL aggregates are functions that run on the server. They can perform statistics such as **Count, Average, Standard Deviation, Min,** and **Max**. Aggregate functions can also be paired with **GROUP BY** to calculate statistics for different groupings of data.

Django supports aggregate functions as of 1.1. Read more about it here: `Django Aggregates <https://docs.djangoproject.com/en/dev/topics/db/aggregation/>`_ ::

    SELECT COUNT(L1), MIN(L1), MAX(L1), AVG(L1), STDDEV(L1) FROM pgd_core_residue GROUP BY aa;

Aggregate functions increase the speed of calculations because they are run on the data in place. Transferring data between the database server and application server requires significant overhead.

------------------------------
Statistics for Dihedral Angles
------------------------------

.. .. image:: directional_average.png

Dihedral angles require special function for average and standard deviation. The special function takes into account that the angles may wrap around from 180 to -180. Both of these functions work like any other aggregate functions. They have also been wrapped in a custom Django Aggregate function to work with django querysets.

-------
Average
-------
::
    IF(DEGREES(ATAN2(
                    -AVG(SIN(RADIANS(ome)))
                    ,-AVG(COS(RADIANS(ome))))

                ) < 0
                ,DEGREES(ATAN2(

                    -AVG(SIN(RADIANS(ome)))
                    ,-AVG(COS(RADIANS(ome))))

                ) + 180
                ,DEGREES(ATAN2(

                    -AVG(SIN(RADIANS(ome)))
                    ,-AVG(COS(RADIANS(ome))))

                ) - 180

            ) AS ome_avg

This works by converting the value into vectors. It adjusts the angles by +180 or -180 depending on whether it is a positive or negative angle. This shifts the vectors into the same space so that they may be averaged.

------------------
Standard Deviation
------------------
::
    SQRT(
            IF (((ome+360)%360 - avgs.ome_avg) < 180
            ,SUM(POW((ome+360)%360-avgs.ome_avg, 2))
            ,SUM(POW(360-((ome+360)%360-avgs.ome_avg),2))

        )/(COUNT(ome)-1))
        AS OME_STDDEV

This function is very similar to a normal standard deviation calculation. The only difference is that a dihedral angle can have two deviations, the short and long way around the circle. We always want to use the shortest distance.

^^^^^^^^^^^^^^^^^
Average Selection
^^^^^^^^^^^^^^^^^

This standard deviation aggregate requires that the average be passed in. There are only two ways to match a list of averages to groups, subqueries or case logic. Neither is an ideal solution but case logic is the lesser of two evils. ::

    CASE SS WHEN 'B' THEN 'foo' WHEN 'H' THEN 'bar' END

This results in queries that are very long (text size), but execution time is fast enough.
