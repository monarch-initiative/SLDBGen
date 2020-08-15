from unittest import TestCase
import os.path
from idg2sl.parsers.entrez_parser import EntrezParser
import idg2sl
from idg2sl import Turner2008Parser

class TestTurner(TestCase):
    """
    This class tests parsing for
    Turner NC, et al., A synthetic lethal siRNA screen identifying genes mediating
    sensitivity to a PARP inhibitor. EMBO J. 2008 May 7;27(9):1368-77. PubMed PMID: 18388863
    We read the first ten lines of the file
    """
    def setUp(self) -> None:
        self.inputfile = os.path.join(os.path.dirname(
            __file__), 'data', 'turner_small.tsv')
        parser = Turner2008Parser(fname=self.inputfile)
        self.turner_list = parser.parse()
        self.first_entry = self.turner_list[0]

    def test_count_entries(self):
        """
        There are 10 unique genes in turner_small.
        There are 10 genes altogether.
        Thus, 10 of the entries should be marked as max
        However, we skip one of the genes (IMPK)
        CDC2 is not an up to date gene symbol and since it is not synthetic lethal,
        the parser skips it. Therefore, we expect 8
        :return:
        """
        self.assertEqual(8, len(self.turner_list))
        num_max = sum([1 for item in self.turner_list if item.is_maximum()])
        self.assertEqual(8, num_max)
