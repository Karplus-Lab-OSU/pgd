

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
        self.indexed_angles = [prefix % f for f in angles]
        self.indexed_fields = [prefix % f for f in fields]
        self.aa_field = prefix % 'aa'
        self.queryset = queryset.values(*(self.indexed_angles+self.indexed_fields+[self.aa_field]))

        self.results = None


    def as_sql(self):
        outer_parts = []
        inner_parts = []
        for field in self.indexed_angles:
            inner_parts.append(
                'ROUND(IF(DEGREES(ATAN2(-AVG(SIN(RADIANS(%(field)s))),-AVG(COS(RADIANS(%(field)s))))) < 0,DEGREES(ATAN2(-AVG(SIN(RADIANS(%(field)s))),-AVG(COS(RADIANS(%(field)s))))) + 180,DEGREES(ATAN2(-AVG(SIN(RADIANS(%(field)s))),-AVG(COS(RADIANS(%(field)s))))) - 180),1) AS avg_%(field)s' % {'field':field}
            )

            outer_parts.append('MIN(%s)' % field)
            outer_parts.append('MAX(%s)' % field)
            outer_parts.append('avg_%s' % field)
            outer_parts.append(
                'SQRT(IF (((%(field)s+360)%%%%360 - avgs.avg_%(field)s) < 180,SUM(POW((%(field)s+360)%%%%360-avgs.avg_%(field)s, 2)),SUM(POW(360-((%(field)s+360)%%%%360-avgs.avg_%(field)s),2)))/(COUNT(%(field)s)-1))AS stddev_%(field)s' % {'field':field}
            )

        for f in self.indexed_fields:
            outer_parts.append('MIN(%(f)s) AS min_%(f)s' % {'f':f})
            outer_parts.append('MAX(%(f)s) as max_%(f)s' % {'f':f})
            outer_parts.append('AVG(%(f)s) as avg_%(f)s' % {'f':f})
            outer_parts.append('STDDEV(%(f)s) as stddev_%(f)s' % {'f':f})

        inner = ','.join(inner_parts)
        outer = ','.join(outer_parts)

        params = {
            'base': self.queryset.query.__str__(),
            'aa_field': self.aa_field,
            'inner': inner,
            'outer': outer
        }

        return 'SELECT pgd_search_segment.%(aa_field)s, %(outer)s FROM (%(base)s) AS pgd_search_segment,(SELECT  %(aa_field)s, %(inner)s FROM (%(base)s) AS pgd_search_segment GROUP BY %(aa_field)s) AS avgs WHERE pgd_search_segment.%(aa_field)s = avgs.%(aa_field)s GROUP BY pgd_search_segment.%(aa_field)s' % params


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
            dict_ = {self.aa_field: row[0] if row[0] else 'total'}

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