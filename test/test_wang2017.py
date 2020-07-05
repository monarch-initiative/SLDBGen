from unittest import TestCase
import os.path
from idg2sl import EntrezParser
import idg2sl


class TestWang(TestCase):
    """
    This tests the function that parses the
    Wang et al, Gene Essentiality Profiling Reveals Gene Networks and Synthetic Lethal Interactions With Oncogenic Ras
    PubMed PMID: 28162770
    which corresponds to the function
    parse_wang_2017()
    We have extracted 30 lines from the input file and put it into
    a test file in test/data/wang2017_small.tsv

    """
    def setUp(self) -> None:
        self.inputfile = os.path.join(os.path.dirname(
            __file__), 'data', 'wang2017_small.tsv')
        self.entrez_file = os.path.join(os.path.dirname(
            __file__), 'data', 'Homo_sapiens.gene_info.gz')
        parser = EntrezParser(self.entrez_file)
        self.wang2017_list = idg2sl.parse_wang_2017(self.inputfile, parser.get_mapping())
        self.first_entry = self.wang2017_list[0]

    def test_count_entries(self):
        """
        There are 22 SL genes in wang2017 small.
        There are 30 genes altogether.
        Thus, 30 of the entries should be marked as max
        :return:
        """
        self.assertEqual(30, len(self.wang2017_list))
        num_max = sum([1 for item in self.wang2017_list if item.is_maximum()])
        self.assertEqual(30, num_max)

    def test_get_symbol(self):
        # Note that the Wanf2017 parse upper cases all gene symbols
        self.assertEqual("NRAS", self.first_entry.get_gene_A_symbol())
        self.assertEqual("SHOC2", self.first_entry.get_gene_B_symbol())

    def test_get_perturbation(self):
        self.assertEqual("mutation", self.first_entry.get_gene_A_pert())
        self.assertEqual("sgRNA", self.first_entry.get_gene_B_pert())

    def test_get_cellosaurus(self):
        self.assertEqual("Ba/F3", self.first_entry.get_cell_line())
        self.assertEqual("CVCL_0161", self.first_entry.get_cellosaurus_id())

    def test_get_cancer_type(self):
        self.assertEqual("", self.first_entry.get_cancer_type())
        self.assertEqual("", self.first_entry.get_ncit_id())

    def test_get_assay(self):
        self.assertEqual("CRISPR-Cas9 Interference assay", self.first_entry.get_assay())

    def test_get_pmid(self):
        self.assertEqual('PMID:28162770', self.first_entry.get_pmid())
