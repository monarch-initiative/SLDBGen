from unittest import TestCase
import os.path
from idg2sl import Luo2009Parser


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
        self.inputfile = os.path.join(os.path.dirname(__file__), 'data', 'luo_small.tsv')
        parser = Luo2009Parser(fname=self.inputfile)
        self.luo_list = parser.parse()
        self.first_entry = self.luo_list[0]

    def test_count_entries(self):
        """
        There are 7 unique genes in Luo small.
        There are 10 genes altogether.
        Thus, 7 of the entries should be marked as max
        :return:
        """
        self.assertEqual(10, len(self.luo_list))
        num_max = sum([1 for item in self.luo_list if item.is_maximum()])
        self.assertEqual(7, num_max)

    def test_get_symbol(self):
        self.assertEqual("KRAS", self.first_entry.get_gene_A_symbol())
        self.assertEqual("ANAPC1", self.first_entry.get_gene_B_symbol())

    def test_get_perturbation(self):
        self.assertEqual("activating_mutation", self.first_entry.get_gene_A_pert())
        self.assertEqual("shRNA", self.first_entry.get_gene_B_pert())

    def test_get_cellosaurus(self):
        self.assertEqual("DLD-1", self.first_entry.get_cell_line())
        self.assertEqual("CVCL_0248", self.first_entry.get_cellosaurus_id())

    def test_get_cancer_type(self):
        self.assertEqual("Colorectal Carcinoma", self.first_entry.get_cancer_type())
        self.assertEqual("NCIT:C2955", self.first_entry.get_ncit_id())

    def test_get_assay(self):
        self.assertEqual("competitive hybridization;multicolor competition assay", self.first_entry.get_assay())

    def test_get_pmid(self):
        self.assertEqual('19490893', self.first_entry.get_pmid())
