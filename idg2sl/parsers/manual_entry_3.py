from idg2sl.parsers.sl_constants import SlConstants
from idg2sl import SyntheticLethalInteraction
from .manual_entry import ManualEntry


class ManualEntry3(ManualEntry):
    """
       In a number of papers, only one or a handful of synthetic lethal interactions are described.
       These are very valuable. It is easiest to enter this information by hand.
       This class is basically the same as ManualEntry but split up to limit the size of the
       entries
       """
    def __init__(self, entrez, ensembl, synonym):
        super().__init__(fname=None, pmid=None, entrez=entrez, ensembl=ensembl, synonym=synonym)
        self.entries = []
        self._add_li_2020()

    def _add_li_2020(self):
        pmid = '32913191'
        braf = 'BRAF'
        cyp2s1 = 'CYP2S1'
        self.create_and_add_sli(geneA=braf, geneB=cyp2s1,
                                geneApert=SlConstants.ACTIVATING_MUTATION, geneBpert=SlConstants.SI_RNA,
                                cell=SlConstants.BCAP_CELL, cellosaurus=SlConstants.BCAP_CELLOSUARUS,
                                ncit=SlConstants.POORLY_DIFFERENTIATED_THYROID_GLAND_CARCINOMA_NCIT,
                                cancer=SlConstants.POORLY_DIFFERENTIATED_THYROID_GLAND_CARCINOMA,
                                assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)



    def get_entries(self):
        return self.entries