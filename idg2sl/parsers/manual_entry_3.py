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
        self.add_bartz_2006()
        self.add_wang2004()
        self.add_rottmann2005()
        self.add_thomas2007()
        self.add_scholl2009()
        self.add_singh2012()
        self.addBajrami()
        self.add_bailey2015()



    def add_bailey2015(self):
        """
        Bailey ML, . Dependence of Human Colorectal Cells Lacking the FBW7 Tumor Suppressor on the Spindle Assembly Checkpoint.
        Genetics. 2015 Nov;201(3):885-95. doi: 10.1534/genetics.115.180653. Epub 2015 Sep 8. PMID: 26354767; PMCID: PMC4649658.
        """
        pmid = '26354767'
        FBXW7 = 'FBXW7'
        # BUBR1 is BUB1B, HGNC:1149
        # MPS1 is IDUA, HGNC:5391
        slgenes = ('BUB1B', 'BUB1', 'IDUA')
        for geneB in slgenes:
            self.create_and_add_sli(geneA=FBXW7, geneB=geneB,
                                    geneApert=SlConstants.KNOCKOUT, geneBpert=SlConstants.RNA_INTERFERENCE_ASSAY,
                                    cell=SlConstants.HCT_116,
                                    cellosaurus=SlConstants.HCT_116_CELLOSAURUS,
                                    assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)

    def addBajrami(self):
        """
        Synthetic lethality of PARP and NAMPT inhibition in triple-negative breast cancer cells.
         EMBO Mol Med. 2012 Oct;4(10):1087-96.  PMID: 22933245; PMC
        """
        pmid = '22933245'
        parp1 = 'PARP1'
        nampt = 'NAMPT'
        self.create_and_add_sli(geneA=parp1, geneB=nampt,
                                geneApert=SlConstants.PHARMACEUTICAL, geneBpert=SlConstants.SI_RNA,
                                assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)

    def add_singh2012(self):
        """
        Singh A, TAK1 inhibition promotes apoptosis in KRAS-dependent colon cancers.
        Cell. 2012 Feb 17;148(4):639-50. doi: 10.1016/j.cell.2011.12.033. PMID: 22341439;
        """
        pmid = '19490892'
        kras = 'KRAS'
        map3k7 = 'MAP3K7'

        self.create_and_add_sli(geneA=kras, geneB=map3k7,
                                geneApert=SlConstants.ACTIVATING_MUTATION, geneBpert=SlConstants.SI_RNA,
                                cell=SlConstants.SW620_CELL, cellosaurus=SlConstants.SW620_CELLOSAURUS,
                                assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)

    def add_scholl2009(self):
        """
        Scholl C, Synthetic lethal interaction between oncogenic KRAS dependency and STK33 suppression in human cancer cells.
        Cell. 2009 May 29;137(5):821-34. doi: 10.1016/j.cell.2009.03.017. PMID: 19490892.
        """
        pmid = '19490892'
        kras = 'KRAS'
        stk33 = 'STK33'
        self.create_and_add_sli(geneA=kras, geneB=stk33,
                                geneApert=SlConstants.ACTIVATING_MUTATION, geneBpert=SlConstants.SI_RNA,
                                assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)

    def add_thomas2007(self):
        """
        Thomas GV,. Hypoxia-inducible factor
        determines sensitivity to inhibitors of mTOR in kidney cancer.
        Nat Med. 2006 Jan;12(1):122-7. doi: 10.1038/nm1337. Epub 2005 Dec 11. PMID: 16341243.
        """
        pmid = '16341243'
        vhl = 'VHL'
        mTOR = 'MTOR'
        self.create_and_add_sli(geneA=vhl, geneB=mTOR,
                                geneApert=SlConstants.SI_RNA, geneBpert=SlConstants.PHARMACEUTICAL,
                                assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)

    def add_rottmann2005(self):
        """
        Rottmann S,  A TRAIL receptor-dependent synthetic lethal relationship between MYC activation and GSK3beta/FBW7 loss of function.
        Proc Natl Acad Sci U S A. 2005 Oct 18;102(42):15195-200. doi: 10.1073/pnas.0505114102.
        PMID: 16210249; PMCID: PMC1257707.
        """
        pmid = '16210249'
        myc = 'MYC'
        gsk3b = 'GSK3B'
        self.create_and_add_sli(geneA=myc, geneB=gsk3b,
                                geneApert=SlConstants.ACTIVATING_MUTATION, geneBpert=SlConstants.SI_RNA,
                                cell=SlConstants.HA1E_CELL, cellosaurus=SlConstants.HA1E_CELLOSAURUS,
                                assay=SlConstants.APOPTOSIS_ASSAY, pmid=pmid)

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

    def add_bartz_2006(self):
        """
        Bartz SR,  Small interfering RNA screens reveal enhanced cisplatin cytotoxicity in tumor cells having both BRCA
        network and TP53 disruptions.
        Mol Cell Biol. 2006 Dec;26(24):9377-86. doi: 10.1128/MCB.01229-06. Epub 2006 Sep 25. PMID: 17000754; PMCID: PMC1698535.
        s shown in Fig. Fig.5A,5A, silencing of BRCA1, BARD1, BRCA2, RAD51, SHFM1, and CHEK1 enhanced growth inhibition
        by cisplatin approximately four- to sevenfold more in the TP53− TOV21G cells than in TP53+ TOV21G cells.
        Note that the current symbol for SHFM1 is SEM1
        """
        pmid = '17000754'
        tp53 = 'TP53'
        sl_genes = ('BRCA1', 'BARD1', 'BRCA2', 'RAD51', 'SEM1', 'CHEK1')
        for geneB in sl_genes:
            self.create_and_add_sli(geneA=tp53, geneB=geneB,
                                    geneApert=SlConstants.KNOCKOUT, geneBpert=SlConstants.SI_RNA,
                                    cell=SlConstants.TOV21G_CELL,
                                    cellosaurus=SlConstants.TOV21G_CELLOSAURUS,
                                    ncit=SlConstants.OVARIAN_CCC_NCIT,
                                    cancer=SlConstants.OVARIAN_CCC,
                                    assay=SlConstants.CYTOTOXICITY_ASSAY, pmid=pmid)

    def add_wang2004(self):
        """
        Wang Y, Engels IH, Knee DA, Nasoff M, Deveraux QL, Quon KC. Synthetic lethal targeting of MYC by activation of
        the DR5 death receptor pathway. Cancer Cell. 2004 May;5(5):501-12. doi: 10.1016/s1535-6108(04)00113-8. PMID: 15144957.
        DR5 is now TNFRSF10B
        """
        pmid = '15144957'
        self.create_and_add_sli(geneA='MYC', geneB='TNFRSF10B',
                                geneApert=SlConstants.ACTIVATING_MUTATION, geneBpert=SlConstants.PHARMACEUTICAL,
                                cell=SlConstants.BJ_CELL,
                                cellosaurus=SlConstants.BJ_CELLOSAURUS,
                                assay=SlConstants.APOPTOSIS_ASSAY, pmid=pmid)

    def get_entries(self):
        return self.entries
