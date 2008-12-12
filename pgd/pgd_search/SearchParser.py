import re

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
    query = Segment.objects
    for res_index,search_res in enumerate(search.residues):
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
            seg_field = 'r'+res_index+'_'+field
            query = query.__dict__['filter' if search_res.__dict__[field+'_include'] else 'exclude'](
                reduce(
                    lambda x,y: x|y,
                    (
                        (
                            Q(**{seg_field+'__gte',float(range_re.split(constraint)[0])}) &
                            Q(**{seg_field+'__lte':float(range_re.split(constraint)[1])})
                        ) if range_re.search(constraint) else (
                            Q(**(
                                {seg_field+"__in" : (aa_choice[1] for aa_index,aa_choice in enumerate(AA_CHOICES) if search_res.aa_int&1<<aa_index)}
                                                                        if field == 'aa_int' else
                                {seg_field        : int(constraint)}    if field == 'terminal_flag' else
                                {seg_field        : float(constraint)}
                            ))
                        ) for constraint in str(segment.__dict__[field]).split(',')
                    )
                ))
    return query
