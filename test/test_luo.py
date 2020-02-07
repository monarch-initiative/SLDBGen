from unittest import TestCase
import os.path
from idg2sl import SyntheticLethalInteraction
from idg2sl import EntrezParser
import idg2sl


class TestLuo(TestCase):
    """
    This tests the function that parses the
    Luo J, et al A genome-wide RNAi screen identifies multiple synthetic lethal
    interactions with the Ras oncogene. Cell. 2009 May 29;137(5):835-48.
    PubMed PMID: 19490893
    which corresponds to the function
    parse_luo2009_supplemental_file_S3
    We have extracted 10 lines from the input file and put it into
    a test file in test/data/luo_small.tsv

    """
    def setUp(self) -> None:
        self.inputfile = os.path.join(os.path.dirname(
            __file__), 'data', 'luo_small.tsv')
        self.entrez_file = os.path.join(os.path.dirname(
            __file__), 'data', 'Homo_sapiens.gene_info.gz')
        parser = EntrezParser(self.entrez_file)
        self.luo_list = idg2sl.parse_luo2009_supplemental_file_S3(self.inputfile, parser.get_mapping())

    def test_entrez(self):
        """
        There are 7 unique genes in Luo small.
        :return:
        """
        self.assertEqual(7, len(self.luo_list))

