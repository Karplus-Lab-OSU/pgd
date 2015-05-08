from django.db.models.aggregates import Aggregate
from django.db.models.sql.aggregates import Aggregate as SQLAggregate

class PGDAggregate(Aggregate):
    """
    Modified to allow Aggregate functions outside of the Django module
    """

class DirectionalStdDev(PGDAggregate):
    alias = 'DirectionalStdDev'
    name =  'DirectionalStdDev'


class DirectionalAvg(PGDAggregate):
    alias = 'DirectionalAvg'
    name =  'DirectionalAvg'


class BinSort(PGDAggregate):
    alias = 'BinSort'
    name =  'BinSort'
