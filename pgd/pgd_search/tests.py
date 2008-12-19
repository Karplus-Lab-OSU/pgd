import unittest
from pgd_search.SearchParser import *
from pgd_search.models import * 
from pgd_core.models import * 
from pgd_splicer.SegmentBuilder import SegmentBuilderTask
from constants import AA_CHOICES, SS_CHOICES
from math import ceil

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
        for i in (1,2,3):

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
    
    
    def testSearchPerSingles(self):
        
        # create Search
        search = Search()
        search.save()
        
        search_residue = Search_residue()
        search_residue.index = 500
        search.residues.add(search_residue)
        search_residue.search = search
        search_residue.a1_include = True
        search_residue.a1 = "1-2"

        # create associated Search_residues (or not)
        for i,j in enumerate(range(int(ceil(1-searchSettings.segmentSize/2.0)),int(ceil(searchSettings.segmentSize/2.0+1.0)))):
            search_residue.index = j
            search_residue.save()

            self.assertEqual(
                # See that the intended query is executed by parse_search
                set(Segment.objects.filter(**{'r%i_a1__gte'%i:1,'r%i_a1__lte'%i:2}).all()),
                set(parse_search(search).all()),
                "Single residue search failed on search index %i"%i
            )

    def testSearchMultiples(self):
        
        # create Search
        search = Search()
        search.save()
        
        search_residue = Search_residue()
        search_residue.index = 500
        search.residues.add(search_residue)
        search_residue.search = search
        search_residue.a1_include = True
        search_residue.a1 = "1-2"

    def testFoo(self):
        ##need at least 1 test or django wont run setUp() this can be removed when we have real tests
        pass

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
            self.assertEqual(validateQueryField(value), False)
