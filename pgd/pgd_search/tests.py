import unittest
from pgd_search.SearchParser import *


class SearchFieldValidationCase(unittest.TestCase):
    def setUp(self):
        pass

    def testSpeaking(self):

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