from unittest import TestCase
import datetime
from selenium import webdriver
from selenium.common.exceptions import NoSuchElementException
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from pgd_search.models import *
from pgd_core.models import *
#from pgd_splicer.SegmentBuilder import SegmentBuilderTask
from pgd_constants import AA_CHOICES, SS_CHOICES
from math import ceil
from search.SearchForm import SearchSyntaxField
import pytz
from django.test import LiveServerTestCase
import sys
import requests
from cStringIO import StringIO

PRO_MIN = -1
PRO_MAX = 3
FIELDS = ['a1','a2','a3','a4','a5','a6','a7','L1','L2','L3','L4','L5','phi','psi','ome','chi1','chi2','chi3','chi4','chi5','bm','bs','bg','occm','occscs','h_bond_energy','zeta']
FIELDS_DICT = {}
for i in range(1, len(FIELDS)+1):
    #shift values into decimels
    FIELDS_DICT[FIELDS[i-1]] = i*.01

class SearchParserValidation(LiveServerTestCase):

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

    def create_protein(self, i, **kwargs):
        protein = Protein()
        protein.code            = '%+i' %i
        protein.threshold       = i
        protein.resolution      = i + .01
        protein.rfactor         = i + .02
        protein.rfree           = i + .03
        protein.pdb_date        = datetime.datetime(2001,1,1, tzinfo=pytz.utc)
        protein.__dict__.update(kwargs)
        protein.save()
        return protein

    def create_chain(self, protein):
        chain = Chain()
        chain.id = '%s%s' % (protein.code, 'A')
        chain.protein = protein
        chain.code = 'A'
        chain.save()
        return chain

    def create_residue(self, i, z, protein, chain, **kwargs):
        residue = Residue()
        residue.protein = protein
        residue.chain = chain
        residue.chainID = chain.code
        residue.chainIndex = z
        #set choice fields.
        residue.aa = self.calculateAA(i)
        residue.ss = self.calculateSS(i)
        for field in FIELDS:
            residue.__dict__[field] = self.calculateResidueField(i, z, field)
        residue.__dict__.update(kwargs)
        residue.save()
        return residue

    def create_bulk_residues(self, i, numResidue, **kwargs):
        protein = self.create_protein(i, **kwargs)
        chain = self.create_chain(protein)
        chainList = [self.create_residue(i, 1, protein, chain) for z in range(numResidue)]
        return chainList


    #creates a set objects with predictable values so we can predict search results
    def setUp(self):
        self.tearDown()

    def tearDown(self):
        Protein.objects.all().delete()

    def testSearchQueryStrings(self):

        # create protein
        protein = self.create_protein(1)
        #create chain
        chain = self.create_chain(protein)
        #create residues
        chainList = [self.create_residue(1, z, protein, chain, a1=z) for z in range(6)]

        # create Search
        search = Search(segmentLength=5)

        # create associated Search_residues
        data = {}
        #data['index'] = 5
        data['residues'] = 1
        data['a1_i_0'] = 1
        data['a1_0'] = "1-5"

        #Set search data equal to search residue parameters
        search.data = data
        search.save()

        self.assertEqual(
            # See that the intended query is executed by parse_search; upper bound of the expected range is NON-INCLUSIVE
            set(chainList[1:6]),
            set(Search.parse_search(search)),
            "Query strings search test failed on range '%s'.  Expected value was %s    Actual values were %s"%(data['a1_0'], set(chainList[1:5]), set(Search.parse_search(search)))
        )

        data['a1_0'] = "1-3"
        search.data = data
        search.save()
        self.assertEqual(
            # See that the intended query is executed by parse_search; upper bound of the expected range is NON-INCLUSIVE
            set(chainList[1:4]),
            set(Search.parse_search(search)),
            "Query strings search test failed on range '%s'"%data['a1_0']
        )


        data['a1_0'] = '<6'
        search.data = data
        self.assertEqual(
            # See that the intended query is executed by parse_search
            set(chainList[0:6]),
            set(Search.parse_search(search)),
            "Query strings search test failed on range '%s'.  Expected value was %s    Actual values were %s"%(data['a1_0'], set(chainList[0:6]), set(Search.parse_search(search)))
        )

        data['a1_0'] = '<=6'
        search.data = data
        self.assertEqual(
            # See that the intended query is executed by parse_search
            set(chainList[0:7]),
            set(Search.parse_search(search)),
            "Query strings search test failed on range '%s'.  Expected value was %s    Actual values were %s"%(data['a1_0'], set(chainList[0:7]), set(Search.parse_search(search)))
        )

        data['a1_0'] = '>5'
        search.data = data
        self.assertEqual(
            # See that the intended query is executed by parse_search
            set(chainList[6:7]),
            set(Search.parse_search(search)),
            "Query strings search test failed on range '%s'.  Expected value was %s    Actual values were %s"%(data['a1_0'], set(chainList[6:7]), set(Search.parse_search(search)))
        )

        data['a1_0'] = '>=5'
        search.data = data
        self.assertEqual(
            # See that the intended query is executed by parse_search
            set(chainList[5:7]),
            set(Search.parse_search(search)),
            "Query strings search test failed on range '%s'.  Expected value was %s    Actual values were %s"%(data['a1_0'], set(chainList[5:7]), set(Search.parse_search(search)))
        )


    def testSearchResolution(self):

        # create protein
        protein1 = self.create_protein(1, resolution=1.1)
        #create chain
        chain = self.create_chain(protein1)
        #create residues
        chainList = [self.create_residue(1, z, protein1, chain) for z in range(3)]

        # create protein
        protein2 = self.create_protein(2, resolution=0.2)
        #create chain
        chain2 = self.create_chain(protein2)
        #create residues
        chainList2 = [self.create_residue(2, i, protein2, chain2) for i in range(3)]

        #create Search and search constraints
        search = Search(segmentLength=0)
        data = {}
        data['residues'] = 0
        data['resolutionMin'] = .1
        data['resolutionMax'] = 1.2

        #Set search data equal to search residue parameters
        search.data = data
        #search.save()
        self.assertEqual(
            # See that both proteins are returned
            set(chainList+chainList2),
            set(Search.parse_search(search)),
            "Resolution search failed on %f-%f   Expected results were %s   Returned results were %s"%(data['resolutionMin'],data['resolutionMax'],list(x.protein.resolution for x in chainList+chainList2),list(y.protein.resolution for y in Search.parse_search(search)))
        )

        data['resolutionMin'] = 0.1
        data['resolutionMax'] = 0.3
        #Set search data equal to search residue parameters
        search.data = data
        #search.save()
        self.assertEqual(
            # See that one protein is returned
            set(chainList2),
            set(Search.parse_search(search)),
            "Resolution search failed on %f-%f   Expected results were %s   Returned results were %s"%(data['resolutionMin'],data['resolutionMax'],list(x.protein.resolution for x in chainList2),list(y.protein.resolution for y in Search.parse_search(search)))
        )

        data['resolutionMin'] = 1
        data['resolutionMax'] = 1.2
        #Set search data equal to search residue parameters
        search.data = data
        search.save()
        self.assertEqual(
            # See that one protein is returned
            set(chainList),
            set(Search.parse_search(search)),
            "Resolution search failed on %f-%f   Expected results were %s   Returned results were %s"%(data['resolutionMin'],data['resolutionMax'],list(x.protein.resolution for x in chainList),list(y.protein.resolution for y in Search.parse_search(search)))
        )

        data['resolutionMin'] = 1.1
        data['resolutionMax'] = 1.1
        #Set search data to exactly the second protein's resolution
        search.data = data
        search.save()
        self.assertEqual(
            # See that one protein is returned
            set(chainList),
            set(Search.parse_search(search)),
            "Resolution search failed on %f-%f   Expected results were %s   Returned results were %s"%(data['resolutionMin'],data['resolutionMax'],list(x.protein.resolution for x in chainList),list(y.protein.resolution for y in Search.parse_search(search)))
        )

        data['resolutionMin'] = .3
        data['resolutionMax'] = .1
        #Set search data to exactly the second protein's resolution
        search.data = data
        search.save()
        self.assertEqual(
            # See that no proteins are returned
            set(),
            set(Search.parse_search(search)),
            "Resolution search failed on %f-%f   Expected results were %s   Returned results were %s"%(data['resolutionMin'],data['resolutionMax'],set(),list(y.protein.resolution for y in Search.parse_search(search)))
        )


    def testSearchThreshold(self):

        #locals = locals()
        chainList = {}
        for t in zip(range(1,6), (1,10,25,25,70)):
            i,thresh = t
            chainList[i] = self.create_bulk_residues( i, 1, threshold = thresh)

        #setattr(obj, name, value)

        # create Search
        search = Search(segmentLength=0)
        data = {}
        data['residues'] = 0
        data['threshold'] = 0
        search.data = data

        self.assertEqual(
            # See that the intended query is executed by parse_search
            set(),
            set(Search.parse_search(search)),
            "Threshold search failed on %i  Expected results were the %s  Returned results were %s"%(data['threshold'],set(), list(y.protein.threshold for y in Search.parse_search(search)))
        )

        data['threshold'] = 25
        search.data = data
        self.assertEqual(
            # See that the intended query is executed by parse_search
            set(chainList[1]+chainList[2]+chainList[3]+chainList[4]),
            set(Search.parse_search(search)),
            "Threshold search failed on %i  Expected results were the %s  Returned results were %s"%(data['threshold'],list(y.protein.threshold for y in chainList[1]+chainList[2]+chainList[3]+chainList[4]), list(y.protein.threshold for y in Search.parse_search(search)))
        )

        data['threshold'] = 75

        self.assertEqual(
            # See that the intended query is executed by parse_search
            set(chainList[1]+chainList[2]+chainList[3]+chainList[4]+chainList[5]),
            set(Search.parse_search(search)),
            "Threshold search failed on %i  Expected results were the %s  Returned results were %s"%(data['threshold'],list(y.protein.threshold for y in chainList[1]+chainList[2]+chainList[3]+chainList[4]+chainList[5]), list(y.protein.threshold for y in Search.parse_search(search)))
        )


        #for index in range(PRO_MIN,PRO_MAX):
            #search.threshold = index
            #search.save()

            #self.assertEqual(
                ## See that the intended query is executed by parse_search
                #set(Segment.objects.filter(protein__threshold=index).all()),
                #set(Search.parse_search(search).all()),
                #"Threshold search failed on %i"%(index)
            #)

    def testSearchCode(self):
        data = {}
        #create residues
        chainList = self.create_bulk_residues( 1, 3, code='gmez')
        chainList2 = self.create_bulk_residues( 2, 3, code='kats')
        # create Search
        search = Search(segmentLength=0)
        search.codes_include = True
        data['proteins_i'] = True
        data['residues'] = 0
        data['proteins'] = 'gmez'
        search.data = data

        self.assertEqual(
            set(chainList),
            set(Search.parse_search(search)),
            "Resolution search failed on %s   Expected results were %s   Returned results were %s"%(data['proteins'],list(x.protein.code for x in chainList),list(y.protein.code for y in Search.parse_search(search)))
        )

        #test inverse
        data['proteins_i'] = False
        self.assertEqual(
            set(chainList2),
            set(Search.parse_search(search)),
            "Resolution search failed on %s   Expected results were %s   Returned results were %s"%(data['proteins'],list(x.protein.code for x in chainList),list(y.protein.code for y in Search.parse_search(search)))
        )

        #testing multi-code and parsing robustness
        data['proteins'] = 'gmez,             kats'
        data['proteins_i'] = True
        self.assertEqual(
            set(chainList+chainList2),
            set(Search.parse_search(search)),
            "Resolution search failed on %s   Expected results were %s   Returned results were %s"%(data['proteins'],list(x.protein.code for x in chainList+chainList2),list(y.protein.code for y in Search.parse_search(search)))
        )

        #test inverse
        data['proteins_i'] = True
        data['proteins'] = 'none'
        self.assertEqual(
            set(),
            set(Search.parse_search(search)),
            "Resolution search failed on the empty set   Expected results were the empty set   Returned results were %s"%(list(y.protein.code for y in Search.parse_search(search)))
        )

        #test inverse
        data['proteins_i'] = False
        data['proteins'] = 'gmez,kats'
        self.assertEqual(
            set(),
            set(Search.parse_search(search)),
            "Resolution search failed on the empty set   Expected results were the empty set   Returned results were %s"%(list(y.protein.code for y in Search.parse_search(search)))
        )

    def testSearchAa(self):
        aa_range = ('a', 'r', 'n', 'd', 'c')
        chainList = {}
        for t in zip(range(1,6), aa_range):
            i,aa_choice = t
            chainList[i] = self.create_bulk_residues( i, 1, aa = aa_choice)

        # create Search
        search = Search(segmentLength=0)
        #search.save()
        data = {}
        data['residues'] = 1

        for aa_index,aa_choice in enumerate(aa_range):
            data['aa_i_0'] = 1
            data['aa_0'] = str(aa_choice)
            search.data = data

            self.assertEqual(
                # See that the intended query is executed by parse_search
                list(chainList[aa_index+1]),
                list(Search.parse_search(search)),
                "Specific AA search failed on '%s'  Expected results were %s  Returned results were %s with %s entries"%(aa_choice, (list(getattr(x, 'aa') for x in (chainList[aa_index+1]))), (list(str(getattr(y, 'aa')) for y in (Search.parse_search(search)))), len(Search.parse_search(search)))
            )

            #test the negation of a single index.
            data['aa_i_0'] = 0
            negatedList = []
            for k in range(1,6):
                if k != (aa_index+1):
                    negatedList += chainList[k]
            search.data = data
            self.assertEqual(
                list(negatedList),
                list(Search.parse_search(search)),
                "Negated AA search failed on '%s'  Expected results were %s  Returned results were %s with %s entries"%(aa_choice, (list(getattr(x, 'aa') for x in (negatedList))), (list(str(getattr(y, 'aa')) for y in (Search.parse_search(search)))), len(Search.parse_search(search)))
            )

        data['aa_i_0'] = 1
        aa_choice = ('a', 'r', 'n')
        data['aa_0'] = (aa_choice)
        search.data = data
        #data['aa_i'] = 0
        self.assertEqual(
            # See that the intended query is executed by parse_search
            list(chainList[1] + chainList[2] + chainList[3]),
            list(Search.parse_search(search)),
            "Specific AA search failed on '%s'  Expected results were %s  Returned results were %s with %s entries"%(aa_choice, (list(getattr(x, 'aa') for x in (chainList[1] + chainList[2] + chainList[3]))), (list(str(getattr(y, 'aa')) for y in (Search.parse_search(search)))), len(Search.parse_search(search)))
        )

        data['aa_i_0'] = 0
        aa_choice = ('a', 'r', 'n')
        data['aa_0'] = (aa_choice)
        search.data = data
        #data['aa_i'] = 0
        self.assertEqual(
            # See that the intended query is executed by parse_search
            list(chainList[4] + chainList[5]),
            list(Search.parse_search(search)),
            "Specific AA search failed on '%s'  Expected results were %s  Returned results were %s with %s entries"%(aa_choice, (list(getattr(x, 'aa') for x in (chainList[4] + chainList[5]))), (list(str(getattr(y, 'aa')) for y in (Search.parse_search(search)))), len(Search.parse_search(search)))
        )


    def testSearchSs(self):

        ss_range = ('H', 'G',  'E', 'T', 'S')
        chainList = {}
        for t in zip(range(1,6), ss_range):
            i,ss_choice = t
            chainList[i] = self.create_bulk_residues( i, 1, ss = ss_choice)

        # create Search
        search = Search(segmentLength=0)
        #search.save()
        data = {}
        data['residues'] = 1

        for ss_index,ss_choice in enumerate(ss_range):
            data['ss_i_0'] = 1
            data['ss_0'] = ss_choice
            search.data = data

            self.assertEqual(
                # See that the intended query is executed by parse_search
                list(chainList[ss_index+1]),
                list(Search.parse_search(search)),
                "Specific ss search failed on '%s'  Expected results were %s  Returned results were %s with %s entries"%(ss_choice, (list(getattr(x, 'ss') for x in (chainList[ss_index+1]))), (list(str(getattr(y, 'ss')) for y in (Search.parse_search(search)))), len(Search.parse_search(search)))
            )

            #test the negation of a single index.
            data['ss_i_0'] = 0
            negatedList = []
            for k in range(1,6):
                if k != (ss_index+1):
                    negatedList += chainList[k]
            search.data = data
            self.assertEqual(
                list(negatedList).sort(),
                list(Search.parse_search(search)).sort(),
                "Negated ss search failed on '%s'  Expected results were %s  Returned results were %s with %s entries"%(ss_choice, (list(getattr(x, 'ss') for x in (negatedList))), (list(str(getattr(y, 'ss')) for y in (Search.parse_search(search)))), len(Search.parse_search(search)))
            )

        data['ss_i_0'] = 1
        ss_choice = ('H', 'G',  'E')
        data['ss_0'] = (ss_choice)
        search.data = data
        #data['ss_i'] = 0
        self.assertEqual(
            # See that the intended query is executed by parse_search
            list(chainList[1] + chainList[2] + chainList[3]).sort(),
            list(Search.parse_search(search)).sort(),
            "Specific ss search failed on '%s'  Expected results were %s  Returned results were %s with %s entries"%(ss_choice, (list(getattr(x, 'ss') for x in (chainList[1] + chainList[2] + chainList[3]))), (list(str(getattr(y, 'ss')) for y in (Search.parse_search(search)))), len(Search.parse_search(search)))
        )

        data['ss_i_0'] = 0
        ss_choice = ('H', 'G',  'E')
        data['ss_0'] = (ss_choice)
        search.data = data
        #data['ss_i'] = 0
        self.assertEqual(
            # See that the intended query is executed by parse_search
            list(chainList[4] + chainList[5]).sort(),
            list(Search.parse_search(search)).sort(),
            "Specific ss search failed on '%s'  Expected results were %s  Returned results were %s with %s entries"%(ss_choice, (list(getattr(x, 'ss') for x in (chainList[4] + chainList[5]))), (list(str(getattr(y, 'ss')) for y in (Search.parse_search(search)))), len(Search.parse_search(search)))
        )

    def testSearchMultipleResidues(self):
        # create protein
        protein = self.create_protein(1)
        #create chain
        chain = self.create_chain(protein)
        #create residues
        chainList = [self.create_residue(1, z, protein, chain) for z in range(1,6)]

        # create Search
        search = Search(segmentLength=0)

        for i, field in enumerate(FIELDS):

            # create associated Search_residues
            data = {}
            data['residues'] = 1
            data['%s_i_-2'%field] = 1
            data['%s_-2'%field] = str(self.calculateResidueField(1, 1, field))
            data['%s_i_0'%field] = 1
            data['%s_0'%field] = str(self.calculateResidueField(1, 3, field))
            data['%s_i_2'%field] = 1
            data['%s_2'%field] = str(self.calculateResidueField(1, 5, field))
            search.data = data
            self.assertAlmostEqual(
                (getattr(chainList[2], '%s'%field)),
                (getattr((Search.parse_search(search).all()[0]), '%s'%field)),#this grabs the first object in the returned search, and attempts to grab the attribute in 'field' out of it
                4,
                "Multiple residue search failed on field %s  Expected result was %s  Returned result was %s"%(field, set((getattr(chainList[2], '%s'%field), )), set(getattr(x, '%s'%field) for x in Search.parse_search(search)))
            )

class SearchFieldValidationCase(LiveServerTestCase):
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
        searchField = SearchSyntaxField()

        for value in validFields:
            self.assertNotEqual(searchField.syntaxPattern.match(value), None, "Valid Field Pattern Failed: '%s'" % value)

        for value in invalidFields:
            self.assertNotEqual(searchField.syntaxPattern.match(value), None)

#Selenium tests
class PersistingSearchOptions(LiveServerTestCase):
    def setUp(self):

       # Create a new instance of the Firefox driver
        self.driver = webdriver.PhantomJS() #webdriver.Firefox()

    def tearDown(self):
        self.driver.quit()

    def test_removed_options_persist(self):
        # Load search page
        self.driver.get(self.live_server_url + "/search")

        # Select the box that indicates number of residues
        residues = self.driver.find_element_by_id("id_residues")
        for option in residues.find_elements_by_tag_name('option'):
            if option.text == "4":
                option.click()

        #Composition
        composition = self.driver.find_element_by_id("id_aa_choices_list_col_2")
        comp_option = composition.find_elements_by_tag_name('li')[0]
        comp_option.click()

        #Conformation
        conformation = self.driver.find_element_by_id("id_ss_choices_list_col_2")
        conf_option = conformation.find_elements_by_tag_name('li')[0]
        conf_option.click()

        conformation_phi = self.driver.find_element_by_id("id_phi_2")
        conformation_phi.clear()
        conformation_phi.send_keys("<=-85,>=85")

        conformation_psi = self.driver.find_element_by_id("id_psi_2")
        conformation_psi.clear()
        conformation_psi.send_keys("<=-80,>=80")

        conformation_omega = self.driver.find_element_by_id("id_ome_2")
        conformation_omega.clear()
        conformation_omega.send_keys("<=-75,>=75")

        #Mobility
        self.driver.find_element_by_id("mobility_header").click()
        mobility = self.driver.find_element_by_id("id_bm_2")
        mobility.clear()
        mobility.send_keys("<35")

        #Angles
        self.driver.find_element_by_id("angles_header").click()
        angles = self.driver.find_element_by_id("id_a1_2")
        angles.send_keys("30")

        #Lengths
        self.driver.find_element_by_id("lengths_header").click()
        lengths = self.driver.find_element_by_id("id_L1_2")
        lengths.send_keys("25")

        #XAngles
        self.driver.find_element_by_id("chi_header").click()
        xangles = self.driver.find_element_by_id("id_chi1_2")
        xangles.send_keys("20")

        #Hackish way to do it, but there doesn't seem to be any other
        #common ways to do it.
        residues = self.driver.find_element_by_id("id_residues")
        for option in residues.find_elements_by_tag_name('option'):
            if option.text == "3":
                option.click()

        for option in residues.find_elements_by_tag_name('option'):
            if option.text == "4":
                option.click()

        #Composition
        self.assertEquals(comp_option.get_attribute("class"), " ")

        #Conformation
        self.assertEquals(conf_option.get_attribute("class"), " ")
        self.assertEquals(conformation_phi.get_attribute("value"), "")
        self.assertEquals(conformation_psi.get_attribute("value"), "")
        self.assertEquals(conformation_omega.get_attribute("value"), "<=-90,>=90")

        #Mobility
        self.assertEquals(mobility.get_attribute("value"), "<25")

        #Angles
        self.assertEquals(angles.get_attribute("value"), "")

        #Lengths
        self.assertEquals(lengths.get_attribute("value"), "")

        #XAngles
        self.assertEquals(xangles.get_attribute("value"), "")

    def test_sidechain_angles_reset(self):
        # Open a new search.
        self.driver.get(self.live_server_url + "/search")

        # Select the amino acid "His" on composition.
        AAs = self.driver.find_element_by_id("id_aa_choices_list_col_0")
        # JMT: There is probably a better way to do this
        his_option = AAs.find_elements_by_tag_name('li')[8]
        his_option.click()

        # Open the "Sidechain Angles" section.
        self.driver.find_element_by_id("sidechain_angle_header").click()

        # Add the filter "<123" to the CbCgCd2 box.
        CbCgCd2_box = self.driver.find_element_by_id("id_HIS__CB_CG_CD2_0")
        CbCgCd2_box.clear()
        CbCgCd2_box.send_keys("<123")

        # Deselect the amino acid "His" on composition.
        his_option.click()

        # Select the amino acid "His" on composition.
        his_option.click()

        # What I expected to see:
        # The CbCgCd2 box should be empty.
        self.assertEquals(CbCgCd2_box.get_attribute("value"), "")

        # What I did see:
        # The CbCgCd2 box contains "<123".
        pass

    def test_sidechain_lengths_reset(self):
        # Open a new search.
        self.driver.get(self.live_server_url + "/search")

        # Select the amino acid "His" on composition.
        AAs = self.driver.find_element_by_id("id_aa_choices_list_col_0")
        # JMT: There is probably a better way to do this
        his_option = AAs.find_elements_by_tag_name('li')[8]
        his_option.click()

        # Open the "Sidechain Lengths" section.
        self.driver.find_element_by_id("sidechain_length_header").click()

        # Add the filter "1" to the CbCg box.
        CbCg_box = self.driver.find_element_by_id("id_HIS__CB_CG_0")
        CbCg_box.clear()
        CbCg_box.send_keys("1")

        # Deselect the amino acid "His" on composition.
        his_option.click()

        # Select the amino acid "His" on composition.
        his_option.click()

        # What I expected to see:
        # The CbCg box should be empty.
        self.assertEquals(CbCg_box.get_attribute("value"), "")

        # What I did see:
        # The CbCg box contains "1".
        pass


class SeleniumTests(LiveServerTestCase):
    fixtures = ['pgd_splicer']

    @classmethod
    def setUpClass(cls):
        cls.driver = webdriver.PhantomJS()
        cls.driver.wait = WebDriverWait(cls.driver, 10)
        super(SeleniumTests, cls).setUpClass()

    @classmethod
    def tearDownClass(cls):
        cls.driver.quit()
        super(SeleniumTests, cls).tearDownClass()

    def test_sidechain_statistics_present(self):
        # Load search page.
        self.driver.get(self.live_server_url + "/search")

        # Run default search.
        self.driver.find_element_by_css_selector('input.submit').click()

        # Confirm that it returns 85 results.
        # JMT: this could be a separate test.
        # Running searches is expensive, though.
        results = self.driver.find_element_by_css_selector('span.results')
        # JMT: if this ever fails, set up WebDriverWait on ajax-loader.png
        self.assertTrue("Results: 85", results.text)

        # Can't save the plot to a file because PhantomJS does not yet
        # support file downloads.

        # Visit statistics link.
        try:
            stats_link = self.driver.wait.until(
                EC.presence_of_element_located((By.LINK_TEXT, 'Statistics'))
            )
            stats_link.click()
        except:
            self.driver.save_screenshot("no-statistics.png")
            self.fail("no stats link")

        # Wait for qtip to disappear.
        # This is problematic!
        qtip_xpath = "//div[contains(., 'Calculating Statistics')]"
        qtip_element = self.driver.find_element_by_xpath(qtip_xpath)
        if qtip_element.is_displayed():
            try:
                qtip_element = self.driver.wait.until(
                    EC.invisibility_of_element_located((By.XPATH, qtip_xpath))
                )
            except:
                self.driver.save_screenshot("no-qtip.png")
                # self.fail("qtip does not disappear")

        # Sidechain statistics are stored in divs named after the residue.
        # Inside each div is two tables -- one for lengths, one for angles.
        # Each table has two rows -- one for mean, one for standard deviation.
        # Each row has one element for each angle or length.

        ss_test_values = {
            "aa_r": {
                "CB_CG": {
                    "avg": "1.523",
                    "stddev": "0.013",
                },
                "CA_CB_CG": {
                    "avg": "112.5",
                    "stddev": "2.4",
                },
            },
        }

        ss_xpath_fmt = "//div[@id='%s']" + \
                       "/table['%d']" + \
                       "/tbody/tr[@class='%s']" + \
                       "/td[@class='%s']"
        ss_h2_fmt = "//div[@id='%s']/h2"

        for residue, elements in ss_test_values.iteritems():
            h2_xpath = ss_h2_fmt % (residue)
            h2_element = self.driver.find_element_by_xpath(h2_xpath)
            for element, rows in elements.iteritems():
                for row, value in rows.iteritems():
                    # Lengths have one _, angles have two.
                    # Lengths are in the first table, angles are in the second.
                    table = row.count('_')
                    ss_xpath = ss_xpath_fmt % (residue, table, row, element)
                    ss_element = self.driver.find_element_by_xpath(ss_xpath)

                    # Value should not be displayed.
                    self.assertFalse(ss_element.is_displayed())

                    # Trigger sidechain.
                    self.driver.find_element_by_xpath(h2_xpath).click()

                    # Value should be displayed and correct.
                    # The field starts as '--' before the value is calculated.
                    self.assertTrue(ss_element.is_displayed())
                    self.driver.wait.until_not(
                        EC.text_to_be_present_in_element((By.XPATH, ss_xpath),
                                                         '--')
                    )
                    self.assertEqual(value, ss_element.text)

                    # Trigger sidechain.
                    self.driver.find_element_by_xpath(h2_xpath).click()

                    # Value should not be displayed.
                    self.assertFalse(ss_element.is_displayed())

    def test_data_dump(self):
        # Load search page.
        self.driver.get(self.live_server_url + "/search")

        # Run default search.
        self.driver.find_element_by_css_selector('input.submit').click()

        # Run data dump.
        dump_link = self.live_server_url + "/search/dump"
        session = requests.Session()
        cookies = self.driver.get_cookies()

        for cookie in cookies:
            session.cookies.set(cookie['name'], cookie['value'])
        response = session.get(dump_link)

        # It should have more than 7 lines.
        self.assertTrue(response.content.count('\n') > 7)


#selenium test for saving the plot
#Works with selenium 2.45.0
class SaveImageAfterSearch(LiveServerTestCase):

    def setUp(self):
        self.driver = webdriver.PhantomJS()

    def tearDown(self):
        self.driver.quit()

    def test_saving_image(self):

         # Load search page
        self.driver.get(self.live_server_url + "/search")

        # Select the box that indicates number of residues
        residues = self.driver.find_element_by_id("id_residues")
        for option in residues.find_elements_by_tag_name('option'):
            if option.text == "4":
                option.click()

        #Clicking the submit button
        self.driver.find_element_by_xpath("//input[@type='submit']").click()

        #Waiting to click for the page to load
        element = WebDriverWait(self.driver, 60).until(
            EC.presence_of_element_located((By.ID, "button-save")))

        #Click the button on the second page
        response = self.driver.find_element_by_id("button-save").click()


class ViewTest(TestCase):
    def home_page_noerror(self):
        response = self.client.get(reverse('/'))
        self.assertEqual(response.status_code, 200)


class CheckDumpTest(TestCase) :


    fixtures = ['pgd_search.json']

    def test_download_tsv(self):

        search = Search(segmentLength=5)

        # create associated Search_residues
        data = {}
        #data['index'] = 5
        data['residues'] = 5
        data['occm_-1'] = 0.58
        data['occm_0'] = 0.58
        data['occm_1'] = 0.58
        data['occscs_-1'] = 0.58
        data['occscs_0'] = 0.58
        data['occscs_1'] = 0.58

        #Set search data equal to search residue parameters
        search.data = data
        search.save()

        from pgd_search.dump.DataDump import Dump
        dump = Dump(search)
        actual = StringIO()
        content_list = []

        try :
            for i in dump:
                actual.write(i)
                #content_list.append(i)
        except IndexError:
            pass
        self.assertIn('0.58\t\t0.58', actual.getvalue())
