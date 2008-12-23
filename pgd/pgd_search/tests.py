import unittest
from pgd_search.models import * 
from pgd_core.models import * 
from pgd_splicer.SegmentBuilder import SegmentBuilderTask
from constants import AA_CHOICES, SS_CHOICES
from math import ceil
from search.views import validateQueryField

PRO_MIN = -1
PRO_MAX = 3
FIELDS = ['a1','a2','a3','a4','a5','a6','a7','L1','L2','L3','L4','L5','phi','psi','ome','chi','bm','bs','bg','h_bond_energy','zeta']
FIELDS_DICT = {}
for i in range(1, len(FIELDS)+1):
    #shift values into decimels
    FIELDS_DICT[FIELDS[i-1]] = i*.01

class SearchParserValidation(unittest.TestCase):
    
    def calculateAA(self, chainIndex):
        #aa_choice = chainIndex-1 if chainIndex-1 < len(AA_CHOICES) else chainIndex-1-len(AA_CHOICES)
        return AA_CHOICES[(chainIndex-1)%len(AA_CHOICES)][0]

    def calculateSS(self, chainIndex):
        #ss_choice = chainIndex-1 if chainIndex-1 < len(SS_CHOICES) else chainIndex-1-len(SS_CHOICES)
        return SS_CHOICES[(chainIndex-1)%len(SS_CHOICES)][0]

    #calculates a value for a field based on other properties
    #   1) whole number indicates protein id
    #   2) first 2 digits after decimal indicate field in FIELDS lookup list/dict
    #   3) last 2 digits indicate chainIndex
    #
    #   example:  1.0305  =>  1.  03   05  =>  protein=1  field=03   chainindex=05
    def calculateResidueField(self, proteinID, chainIndex, field):
        return proteinID + FIELDS_DICT[field] + (chainIndex*.0001)

    #creates a set objects with predictable values so we can predict search results
    def setUp(self):
        # create a series of proteins, chains, segments
        for i in range(PRO_MIN,PRO_MAX):

            # create protein
            protein = Protein()
            protein.code            = '%+i' %i
            protein.threshold       = i
            protein.resolution      = i + .01
            protein.rfactor         = i + .02

            protein.save()

            #create chains
            for chainID in ['A','B']:
                chain = Chain()
                chain.id = '%s%s' % (protein.code, chainID)
                chain.protein = protein
                chain.code = chainID
                chain.save()

                # create 15 residues per chain
                for z in range(1, 16):
                    residue = Residue()
                    residue.protein = protein
                    residue.chain = chain
                    residue.chainID = chainID
                    residue.chainIndex = z

                    #set choice fields.  
                    residue.aa = self.calculateAA(z)
                    residue.ss = self.calculateSS(z)

                    for field in FIELDS:
                        residue.__dict__[field] = self.calculateResidueField(i, z, field)

                    residue.save()

        #core tables are populated, run SegmentBuilder to populate segments
        print
        builder = SegmentBuilderTask('Test Suite Builder')
        builder._work(None)
    
    def testSearchQueryStrings(self):
        
        # create Search
        search = Search(segmentLength=0)
        search.save()
        
        search_residue = Search_residue()
        search_residue.index = -4
        search.residues.add(search_residue)
        search_residue.search = search
        search_residue.a1_include = True
    
        # create associated Search_residues (or not)
        for min,max in [(x,y) for x in range(PRO_MIN,PRO_MAX) for y in range(x+1,PRO_MAX)]:
            search_residue.a1 = "%g-%g"%(min,max)
            search_residue.save()
            print Search.parse_search(search).all().count()
            self.assertEqual(
                # See that the intended query is executed by parse_search
                set(Segment.objects.filter(**{'r0_a1__gte':min,'r0_a1__lte':max}).all()),
                set(Search.parse_search(search).all()),
                "Query strings search test failed on range '%s'"%search_residue.a1
            )
    
    def testSearchResolution(self):
        
        # create Search
        search = Search(segmentLength=0)

        for min,max in [(x,y) for x in range(PRO_MIN,PRO_MAX) for y in range(x+1,PRO_MAX)]:
            search.resolution_min = min
            search.resolution_max = max
            search.save()

            self.assertEqual(
                # See that the intended query is executed by parse_search
                set(Segment.objects.filter(protein__resolution__gte=min,protein__resolution__lte=max).all()),
                set(Search.parse_search(search).all()),
                "Resolution search failed on %i-%i"%(min,max)
            )

    def testSearchThreshold(self):
        
        # create Search
        search = Search(segmentLength=0)

        for index in range(PRO_MIN,PRO_MAX):
            search.threshold = index
            search.save()

            self.assertEqual(
                # See that the intended query is executed by parse_search
                set(Segment.objects.filter(protein__threshold=index).all()),
                set(Search.parse_search(search).all()),
                "Threshold search failed on %i"%(index)
            )

    def testSearchAa(self):
        
        # create Search
        search = Search(segmentLength=0)
        search.save()

        search_residue = Search_residue()
        search_residue.index = -4
        search.residues.add(search_residue)
        search_residue.search = search
        search_residue.aa_int_include = True
        search_residue.aa_int = 0
        for aa_index,aa_choice in enumerate(AA_CHOICES):
            search_residue.aa_int = 1<<aa_index
            search_residue.save()
            
            self.assertEqual(
                # See that the intended query is executed by parse_search
                set(Segment.objects.filter(r0_aa=aa_choice[0]).all()),
                set(Search.parse_search(search).all()),
                "Specific AA search failed on '%s'"%(aa_choice[1])
            )
        
        search_residue.aa_int = 0
        
        for aa_index,aa_choice in enumerate(AA_CHOICES):
            search_residue.aa_int += 1<<aa_index
            search_residue.save()
            self.assertEqual(
                # See that the intended query is executed by parse_search
                set(Segment.objects.filter(r0_aa__in=[aa_c[0] for aa_i,aa_c in enumerate(AA_CHOICES) if search_residue.aa_int&1<<aa_i]).all()),
                set(Search.parse_search(search).all()),
                "Multiple AA search failed on '%s'"%'\', \''.join([aa_c[1] for aa_i,aa_c in enumerate(AA_CHOICES) if search_residue.aa_int&1<<aa_i])
            )
    
    def testSearchMultipleResidues(self):
        
        # create Search
        search = Search(segmentLength=0)
        search.save()

        for i,j in enumerate(range(int(ceil(1-searchSettings.segmentSize/2.0)),int(ceil(searchSettings.segmentSize/2.0+1.0)))):
            search_residue = Search_residue()
            search_residue.index = j
            search.residues.add(search_residue)
            search_residue.search = search
            search_residue.a1_include = True
            search_residue.a1 = "1-2"
            search_residue.save()

            self.assertEqual(
                # See that the intended query is executed by parse_search
                set(Segment.objects.filter(
                    reduce(
                        lambda x,y: x&y,
                        (Q(**{'r%i_a1__gte'%k:1,'r%i_a1__lte'%k:2})
                        for k in range(i+1))
                    )
                ).all()),
                set(Search.parse_search(search).all()),
                "Multiple residue search failed on %i residues"%(i+1)
            )
    
    def testSearchMultipleFields(self):
        
        # create Search
        search = Search(segmentLength=0)
        search.save()
        search_residue = Search_residue()
        search_residue.index = -4
        search.residues.add(search_residue)
        search_residue.search = search

        for prot_dicts,filter_dict,search_dict in    [
                        (
                            True,
                            {'protein__resolution__gte': 1, 'protein__resolution__lte': 2},
                            {'resolution_min': 1, 'resolution_max': 2}
                        ),
                        (
                            True,
                            {'protein__threshold': 1},
                            {'threshold': 1}
                        )
                    ] + [(
                        False,
                        {'r0_' + field + '__gte': 1, 'r0_' + field + '__lte': 2},
                        {field+'_include': True, field: '1-2'}
                    ) for field in ('a1','a2','a3','a4','a5','a6','a7','L1','L2','L3','L4','L5','phi','psi','ome','chi','bm','bs','bg','h_bond_energy','zeta')]:
            (search if prot_dicts else search_residue).__dict__.update(search_dict)

            search.save()                
            search_residue.save()
            self.assertEqual(
                # See that the intended query is executed by parse_search
                set(Segment.objects.filter(**filter_dict).all()),
                set(Search.parse_search(search).all()),
                "Multiple fields search failed",
            )

class SearchFieldValidationCase(unittest.TestCase):
    def setUp(self):
        pass

    def testFieldSyntaxParser(self):
        validFields = [
            '1',
            '1-1',
            '1,2,3',
            '1-2',
            '1-1',
            '1-2,5-6',
            '1,23,456',
            '1-23',
            '1-23,456-7890',
            '-1',
            '-23',
            '.5',
            '.123',
            '0.5',
            '0.5,0.6',
            '0.5-0.6',
            '0.5-0.6,0.5',
            '0.5-123.123',
            '-1-2',
            '1--2',
            '-1--2'
        ]

        invalidFields = []
        
        for value in validFields:
            self.assertEqual(validateQueryField(value), True, "Valid Field Pattern Failed: '%s'" % value)

        for value in invalidFields:
            self.assertEqual(validateQueryField(value), None)
