from unittest import TestCase
import os.path
from idg2sl import Steckel2012Parser
from idg2sl import SlConstants


class TestSteckel2012(TestCase):
    """
    This tests the function that parses the
    Steckel M, et al. Determination of synthetic lethal interactions in KRAS oncogene-dependent cancer cells reveals novel
    therapeutic targeting strategies. Cell Res. 2012 Aug;22(8):1227-45.  PubMed PMID: 22613949
    We have extracted 20 lines from the input file and put it into
    a test file in test/data/luo_small.tsv. THe top 10 are hits and the bottom ten are not

    """
    def setUp(self) -> None:
        self.inputfile = os.path.join(os.path.dirname(
            __file__), 'data', 'steckel-2012-small.tsv')
        parser = Steckel2012Parser(fname=self.inputfile)
        self.steckel_list = parser.parse()
        self.first_entry = self.steckel_list[0]

    def test_count_entries(self):
        """
        There are 20 unique genes in steckel-2012-small.
        There are 20 genes altogether. One of the genes is KRAS itself,
        which we filter out. Therefore, we have 19 for both tests
        :return:
        """
        expected_genes = 19
        self.assertEqual(expected_genes, len(self.steckel_list))
        num_max = sum([1 for item in self.steckel_list if item.is_maximum()])
        self.assertEqual(expected_genes, num_max)

    def test_get_symbol(self):
        self.assertEqual("KRAS", self.first_entry.get_gene_A_symbol())
        self.assertEqual("POLR2A", self.first_entry.get_gene_B_symbol())

    def test_get_perturbation(self):
        self.assertEqual("activating mutation", self.first_entry.get_gene_A_pert())
        self.assertEqual("siRNA", self.first_entry.get_gene_B_pert())

    def test_get_cellosaurus(self):
        self.assertEqual(SlConstants.HCT_116, self.first_entry.get_cell_line())
        self.assertEqual(SlConstants.HCT_116_CELLOSAURUS, self.first_entry.get_cellosaurus_id())
