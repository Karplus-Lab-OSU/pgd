class DirectionalStatisticsQuery():
    """
    This is a specialized query that uses aggregates and subqueries to
    efficiently calculate directional mean and directional standard deviation
    in a single query.
    """

    def __init__(self, angles, fields, prefix, queryset):
        self.angles = angles
        self.fields = fields
        self.combined = angles+fields

        # create set of fields to query.  these require the django style
        # prefix for the fields.
        self.queryset = queryset.values(*([prefix % f for f in angles] + 
                                        [prefix % f for f in fields] + 
                                        [prefix % 'aa']))
        self.results = None

    def as_sql(self):       
        outer_parts = []
        inner_parts = []
        for field in self.angles:
            inner_parts.append(
                'ROUND(IF(DEGREES(ATAN2(-AVG(SIN(RADIANS(%(f)s))),-AVG(COS(RADIANS(%(f)s))))) < 0,DEGREES(ATAN2(-AVG(SIN(RADIANS(%(f)s))),-AVG(COS(RADIANS(%(f)s))))) + 180,DEGREES(ATAN2(-AVG(SIN(RADIANS(%(f)s))),-AVG(COS(RADIANS(%(f)s))))) - 180),1) AS avg_%(f)s' % {'f':field}
            )

            outer_parts.append('MIN(%s)' % field)
            outer_parts.append('MAX(%s)' % field)
            outer_parts.append('avg_%s' % field)
            outer_parts.append(
                'SQRT(IF (((%(f)s+360)%%%%360 - avgs.avg_%(f)s) < 180,SUM(POW((%(f)s+360)%%%%360-avgs.avg_%(f)s, 2)),SUM(POW(360-((%(f)s+360)%%%%360-avgs.avg_%(f)s),2)))/(COUNT(%(f)s)-1))AS stddev_%(f)s' % {'f':field}
            )

        for f in self.fields:
            f = {'f':f}
            outer_parts.append('MIN(%(f)s) AS min_%(f)s' % f)
            outer_parts.append('MAX(%(f)s) as max_%(f)s' % f)
            outer_parts.append('AVG(%(f)s) as avg_%(f)s' % f)
            outer_parts.append('STDDEV(%(f)s) as stddev_%(f)s' % f)

        inner = ','.join(inner_parts)
        outer = ','.join(outer_parts)

        params = {
            'base': self.queryset.query.__str__(),
            'inner': inner,
            'outer': outer
        }

        return 'SELECT residues.aa, %(outer)s FROM (%(base)s) AS residues, \
            (SELECT aa, %(inner)s FROM (%(base)s) AS residues GROUP BY aa) \
            AS avgs WHERE residues.aa = avgs.aa GROUP BY residues.aa \
            WITH ROLLUP' % params

    def __str__(self):
        if not self.results:
            self._execute()

        return self.results.__str__()


    def execute(self):
        if self.results:
            return self.results

        return self._execute()


    def _execute(self):
        """
        Private method for executing query, always runs query and then updates cache
        """
        from django.db import connection, transaction
        cursor = connection.cursor()
        cursor.execute(self.as_sql())
        results = []

        # iterate results and process results into a dictionary
        row = cursor.fetchone()
        combined = self.combined
        while row:
            dict_ = {'aa': row[0] if row[0] else 'total'}

            for i in range(len(combined)):
                f = combined[i]
                dict_['min_%s' % f] = row[4*i+1]
                dict_['max_%s' % f] = row[4*i+2]
                dict_['avg_%s' % f] = row[4*i+3]
                dict_['stddev_%s' % f] = row[4*i+4]

            results.append(dict_)
            row = cursor.fetchone()

        self.results = results
        return results


    def __iter__(self):
        if not self.results:
            self._execute()

        return enumerate(self.results)