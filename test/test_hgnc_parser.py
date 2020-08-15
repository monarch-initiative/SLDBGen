from unittest import TestCase
import os.path
from idg2sl import HgncParser




class TestHgncParser(TestCase):
    """
    This class tests the HgncParser
    """
    def setUp(self) -> None:
        self.inputfile = os.path.join(os.path.dirname(
            __file__), 'data', 'hgnc_small.txt')
        hgncparser = HgncParser(self.inputfile)
        self.d = hgncparser.get_dictionary() 
       
    def test_count_entries(self):
        #hgncparser = HgncParser(self.inputfile)
         self.assertEqual(10,len(self.d))
