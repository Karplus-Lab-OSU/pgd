import re
from pgd_search.models import Segment,searchSettings
from math import ceil
import re
from django.db.models import Q


range_re = re.compile("(?<=[^-])-")

"""
Validates a field to make sure that it has valid syntax for a search field
"""
def validateQueryField(str):
    r = re.compile(r'^(-?([1-9]\d*|0)(\.\d+)?|(\.\d+))(-(-?([1-9]\d*|0)(\.\d+)?|(\.\d+)))?(,(-?([1-9]\d*|0)(\.\d+)?|(\.\d+))(-(-?([1-9]\d*|0)(\.\d+)?|(\.\d+)))?)*$')
    m = r.match(str)
    return not m == None

"""
Parses a search into a Django query
"""
def parse_search(search):
    search_codes = (x.code for x in search.codes.all())
    query = Segment.objects 
    if search.codes_include:
        query = query.__getattribute__('filter' if search.codes_include else 'exclude')(protein__in=search_codes)
    for search_res in search.residues.all():
        for field in filter(
                lambda x: search_res.__dict__[x+'_include'] != None,
                (
                    'aa_int',
                    'a1',   'a2',   'a3',   'a4',   'a5',   'a6',   'a7',
                    'L1',   'L2',   'L3',   'L4',   'L5',
                    'ss',
                    'phi',  'psi',  'ome',  'chi',
                    'bm',   'bs',   'bg',
                    'h_bond_energy',
                    'zeta',
                    'terminal_flag',
                )):
            seg_field = 'r%i_%s'%((search_res.index+int(ceil(searchSettings.segmentSize/2.0)-1)),field)
            query = query.__getattribute__('filter' if search_res.__dict__[field+'_include'] else 'exclude')(
                reduce(
                    lambda x,y: x|y,
                    (
                        (
                            Q(**{seg_field+'__gte' : float(range_re.split(constraint)[0])}) &
                            Q(**{seg_field+'__lte' : float(range_re.split(constraint)[1])})
                        ) if range_re.search(constraint) else (
                            Q(**(
                                {seg_field+"__in"  : (aa_choice[1] for aa_index,aa_choice in enumerate(AA_CHOICES) if search_res.aa_int&1<<aa_index)}
                                                                        if field == 'aa_int' else
                                {seg_field         : int(constraint)}    if field == 'terminal_flag' else
                                {seg_field         : float(constraint)}
                            ))
                        ) for constraint in str(search_res.__dict__[field]).split(',')
                    )
                ))
    return query
