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
        self._add_ali_2020()
        self._add_ertay_2020()
        self._add_neggers_2020()

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

    def _add_ali_2020(self):
        """
        ATR (AZD6738), ATM (AZ31) or Wee1 (AZD1775) monotherapy was selectively toxic in XRCC1 deficient cells.
        Selective synergistic toxicity was evident when olaparib was combined with AZD6738, AZ31 or AZD1775. The most
        potent synergistic interaction was evident with the AZD6738 and olaparib combination therapy. In clinical
        cohorts, ATR, ATM or Wee1 overexpression in XRCC1 deficient breast cancer was associated with poor outcomes.
        """
        pmid = '33425022'
        xrcc1 = 'XRCC1'
        sl_genes = {'ATR', 'ATM', 'WEE1'}
        parp1 = 'PARP1'
        parp1_id = self.get_ncbigene_curie(parp1)
        for geneB in sl_genes:
            self.create_and_add_sli(geneA=xrcc1, geneB=geneB,
                                    geneApert=SlConstants.LOF_MUTATION, geneBpert=SlConstants.PHARMACEUTICAL,
                                    cell='MDA-MB-231_XRCC1_KO',
                                    background_dependency_gene_symbol=parp1,
                                    background_dependency_gene_id=parp1_id,
                                    ncit=SlConstants.BREAST_CARCINOMA_NCIT,
                                    cancer=SlConstants.BREAST_CARCINOMA,
                                    assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)

    def _add_ertay_2020(self):
        pmid = '33221821'
        pten = 'PTEN'
        wdhd1 = 'WDHD1'
        self.create_and_add_sli(geneA=pten, geneB=wdhd1,
                                geneApert=SlConstants.LOF_MUTATION, geneBpert=SlConstants.SI_RNA,
                                cell=SlConstants.MDAMB468_CELL,
                                cellosaurus=SlConstants.MDAMB468_CELLOSAURUS,
                                ncit=SlConstants.TRIPLE_NEG_BREAST_CARCINOMA_NCIT,
                                cancer=SlConstants.TRIPLE_NEG_BREAST_CARCINOMA,
                                assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)

    def _add_neggers_2020(self):
        pmid = '33326793'
        vps4a = 'VPS4A'
        vps4b = 'VPS4B'
        self.create_and_add_sli(geneA=vps4a, geneB=vps4b,
                                geneApert=SlConstants.CRISPR_CAS9, geneBpert=SlConstants.LOF_MUTATION,
                                cell=SlConstants.RMS_CELL,
                                cellosaurus=SlConstants.RMS_CELLOSAURUS,
                                ncit=SlConstants.EMBRYONAL_RHABDOMYOSARCOMA_NCIT,
                                cancer=SlConstants.EMBRYONAL_RHABDOMYOSARCOMA,
                                background_dependency_status='16q loss',
                                assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)

    def get_entries(self):
        return self.entries
