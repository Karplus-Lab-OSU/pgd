from django.db.models.aggregates import Aggregate
from django.db.models.sql.aggregates import Aggregate as SQLAggregate

class PGDAggregate(Aggregate):
    """
    Modified to allow Aggregate functions outside of the Django module
    """

    def add_to_query(self, query, alias, col, source, is_summary):
        """Add the aggregate to the nominated query.

        This method is used to convert the generic Aggregate definition into a
        backend-specific definition.

         * query is the backend-specific query instance to which the aggregate
           is to be added.
         * col is a column reference describing the subject field
           of the aggregate. It can be an alias, or a tuple describing
           a table and column name.
         * source is the underlying field or aggregate definition for
           the column reference. If the aggregate is not an ordinal or
           computed type, this reference is used to determine the coerced
           output type of the aggregate.
         * is_summary is a boolean that is set True if the aggregate is a
           summary value rather than an annotation.
        """
        klass = globals()['%sSQL' % self.name]
        aggregate = klass(col, source=source, is_summary=is_summary, **self.extra)
        
        # Validate that the backend has a fully supported, correct
        # implementation of this aggregate
        query.aggregates[alias] = aggregate
        self.aggregate = aggregate
    

class DirectionalStdDev(PGDAggregate):
    alias = 'DirectionalStdDev'
    name =  'DirectionalStdDev'


class DirectionalStdDevSQL(SQLAggregate):
    is_computed = True
    sql_function = ''
    sql_template =  '%(function)sSQRT(IF((MOD(%(field)s+360,360)-%(avg)s)<180,SUM(POW(MOD(%(field)s+360,360)-%(avg)s,2)),SUM(POW(360-(MOD(%(field)s+360,360)-%(avg)s),2)))/(COUNT(%(field)s)-1))'


class DirectionalAvg(PGDAggregate):
    alias = 'DirectionalAvg'
    name =  'DirectionalAvg'


class DirectionalAvgSQL(SQLAggregate):
    is_computed = True
    sql_function = ''
    sql_template = '%(function)sIF(DEGREES(ATAN2(-AVG(SIN(RADIANS(%(field)s))),-AVG(COS(RADIANS(%(field)s))))) < 0,DEGREES(ATAN2(-AVG(SIN(RADIANS(%(field)s))),-AVG(COS(RADIANS(%(field)s))))) + 180,DEGREES(ATAN2( -AVG(SIN(RADIANS(%(field)s))),-AVG(COS(RADIANS(%(field)s))))) - 180)'


class BinSort(PGDAggregate):
    alias = 'BinSort'
    name =  'BinSort'


class BinSortSQL(SQLAggregate):
    sql_function = ''
    sql_template = '%(function)sFLOOR((IF(%(field)s<%(offset).16f,360,0)+%(field)s-%(offset).16f)/%(bincount).16f)-IF(%(field)s=%(max).16f,1,0)'
