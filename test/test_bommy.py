from unittest import TestCase
import os.path
from idg2sl import Bommi2008Parser
import idg2sl



class TestBommy(TestCase):
    """
    This class tests parsing for
     Bommi-Reddy A, et al  Kinase requirements in human cells: III. Altered kinase
    requirements in VHL-/- cancer cells detected in a pilot synthetic lethal screen.
    Proc Natl Acad Sci U S A. 2008 Oct 28;105(43):16484-9. PubMed PMID: 18948595.
    We read the first ten lines of the file
    """

    def setUp(self) -> None:
        self.inputfile = os.path.join(os.path.dirname(
            __file__), 'data', 'bommy_small.tsv')
        parser = Bommi2008Parser(fname=self.inputfile)
        self.bommy_list = parser.parse()
        self.first_entry = self.bommy_list[0]

    def test_count_entries(self):
        """
        There are 2 unique genes in Luo small.
        There are 10 genes altogether.
        Thus, 2 of the entries should be marked as max
        :return:
        """
        self.assertEqual(10, len(self.bommy_list))
        num_max = sum([1 for item in self.bommy_list if item.is_maximum()])
        self.assertEqual(2, num_max)

    def test_get_cancer_type(self):
        self.assertEqual("Clear Cell Renal Cell Carcinoma", self.first_entry.get_cancer_type())
        self.assertEqual("NCIT:C4033", self.first_entry.get_ncit_id())