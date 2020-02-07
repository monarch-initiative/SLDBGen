import unittest
from idg2sl import SyntheticLethalInteraction


class TestSyntheticLethalInteraction(unittest.TestCase):
    def setUp(self) -> None:
        self.parameters = {
             'gene_A_symbol': 'symA',
             'gene_A_id': 'idA',
             'gene_B_symbol': 'symB',
             'gene_B_id': 'idB',
             'gene_A_pert': 'pert1',
             'gene_B_pert': 'pert2',
             'effect_type': 'thisEff',
             'effect_size': '20',
             'cell_line': 'cellLine42',
             'cellosaurus_id': 'csID1',
             'cancer_type': 'melanoma',
             'ncit_id': 'ncit1234',
             'assay': 'thisAssay',
             'pmid': '27453043',
             'SL': True
        }
        self.no_getter = []
        self.sli = SyntheticLethalInteraction(**self.parameters)

    def test_attributes_and_getters(self):
        for param, value in self.parameters.items():
            if param in self.no_getter:
                continue
            with self.subTest():
                self.assertEqual(getattr(self.sli, param), value,
                                 msg="Attribute {} has wrong value ({} != {})"
                                 .format(param, getattr(self.sli, param), value))

                getter_name = "get_" + param
                self.assertTrue(hasattr(self.sli, getter_name),
                                msg="No getter method for {}!".format(param))
                getter = getattr(self.sli, getter_name)
                self.assertEqual(getter(), value,
                                 msg="Getter for {} returned the wrong value ({} != {})"
                                 .format(param, getter(), value))


