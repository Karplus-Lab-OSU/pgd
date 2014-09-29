from django.db.models import Min, Max, Avg, StdDev

from pgd_search.statistics.aggregates import DirectionalAvg, DirectionalStdDev


class DirectionalStatisticsQuery():
    """
    This is a specialized query that performs two queries and merges the
    results.  The first calculates any aggregate functions that can be completed
    in a normal query.  min/max/avg for all fields, and stddev for linear
    fields.  The second query calculates directional stddev for any fields that
    require it using the avgs calculated in the first query
    """

    def __init__(self, angles, fields, prefix, queryset):
        self.angles = angles
        self.fields = fields
        self.combined = angles+fields
        self.prefix = prefix
        # create set of fields to query.  these require the django style
        # prefix for the fields.
        self.queryset = queryset
        self.results = None

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
        annotations = {}
        aa_rows = {}
        p = self.prefix
        aa_field = p%'aa'
        
        # main aggregate functions
        for field in self.angles:
            annotations['min_%s' % field] = Min(p%field)
            annotations['max_%s' % field] = Max(p%field)
            annotations['avg_%s' % field] = DirectionalAvg(p%field)
        for field in self.fields:
            annotations['min_%s' % field] = Min(p%field)
            annotations['max_%s' % field] = Max(p%field)
            annotations['avg_%s' % field] = Avg(p%field)
            annotations['stddev_%s' % field] = StdDev(p%field, sample=True)
        
        # query with all aggregate values that can be calculated in a standard
        # query.  save query in a list so that its members can be modified
        query = self.queryset
        query = query.values(p%'aa')
        query = query.annotate(**annotations)
        
        results = list(query)
        
        # construction 2nd query for DirectionStdDev calculations for each
        # dihedral angle.  only needed if there is at least one dihedral angle
        if self.angles:
            annotations = {}
            for row in results:
                aa_rows[row[aa_field]] = row
                for field in self.angles:
                    avg=row['avg_%s' % field]
                    if avg:
                        annotations['stddev_%s' % field] = DirectionalStdDev(p%field, avg=avg)
                    else:
                        row['stddev_%s' % field] = None
            
            if annotations:
                query = self.queryset
                query = query.values(aa_field)
                query = query.annotate(**annotations)
            
            # update the original results with the results of the 2nd query
            for row in query:
                outer_row = aa_rows[row[aa_field]]
                outer_row.update(row)
        
        # rename aa columns
        if self.prefix != '%s':
            for row in results:
                row['aa'] = row[aa_field]
                del row[aa_field]
        
        self.results = results
        return results

    def __iter__(self):
        if not self.results:
            self._execute()
        return iter(self.results)


class DirectionalStatisticsTotalQuery(DirectionalStatisticsQuery):
    """
    Extension of class for doing totals
    """
    def _execute(self):
        """
        Private method for executing query, always runs query and then updates cache
        """
        annotations = {}
        aa_rows = {}
        p = self.prefix
        aa_field = p%'aa'
        
        # main aggregate functions
        for field in self.angles:
            annotations['min_%s' % field] = Min(p%field)
            annotations['max_%s' % field] = Max(p%field)
            annotations['avg_%s' % field] = DirectionalAvg(p%field)
        for field in self.fields:
            annotations['min_%s' % field] = Min(p%field)
            annotations['max_%s' % field] = Max(p%field)
            annotations['avg_%s' % field] = Avg(p%field)
            annotations['stddev_%s' % field] = StdDev(p%field, sample=True)
        
        # query with all aggregate values that can be calculated in a standard
        # query.  save query in a list so that its members can be modified
        query = self.queryset
        totals = query.aggregate(**annotations)
        
        # construction 2nd query for DirectionStdDev calculations for each
        # dihedral angle.
        if self.angles:
            annotations = {}
            for field in self.angles:
                avg=totals['avg_%s' % field]
                if avg:
                    annotations['stddev_%s' % field] = DirectionalStdDev(field, avg=avg)
                else:
                    totals['stddev_%s' % field] = None
            
            if annotations:
                query = self.queryset
                stddevs = query.aggregate(**annotations)
                totals.update(stddevs)
        
        totals['aa'] = 'total'
        self.results = [totals]
        
        return self.results
