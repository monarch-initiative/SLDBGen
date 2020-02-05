import unittest
from utils.entrez_lookup import EntrezLookup
from collections import defaultdict


class TestUtils(unittest.TestCase):

    def test_entrez_lookup_type(self):
        el = EntrezLookup()
        self.assertTrue(hasattr(el, 'reverse_lookup'))
        self.assertTrue(isinstance(el.reverse_lookup, defaultdict))

