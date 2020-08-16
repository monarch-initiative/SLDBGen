from idg2sl.parsers.sl_constants import SlConstants
from idg2sl import SyntheticLethalInteraction

breastCarcinoma = "Breast Carcinoma"
ncitBreastCarcinoma = "NCIT:C4872"


class ManualEntry:
    """
    In a number of papers, only one or a handful of synthetic lethal interactions are described.
    These are very valuable. It is easiest to enter this information by hand.
    """

    def __init__(self):
        self.entries = []
        self._add_reid_2016()
        self._add_parthak_2015()
        self._add_sultana_2013()
        self._add_hoffmann_2014()
        self._add_scholl_2009()
        self._add_barbie_2009()
        self._add_chen_2020()
        self._add_yamada_2020()
        self._add_lelij_2020()
        self._add_liu_2020()

    def create_sli(self, geneA, geneAid, geneB, geneBid, geneApert, geneBpert, effecttype, effectsize, cell, cellosaurus,
                   cancer, ncit, assay, pmid):
        sli = SyntheticLethalInteraction(gene_A_symbol=geneA,
                                         gene_A_id=geneAid,
                                         gene_B_symbol=geneB,
                                         gene_B_id=geneBid,
                                         gene_A_pert=geneApert,
                                         gene_B_pert=geneBpert,
                                         effect_type=effecttype,
                                         effect_size=effectsize,
                                         cell_line=cell,
                                         cellosaurus_id=cellosaurus,
                                         cancer_type=cancer,
                                         ncit_id=ncit,
                                         assay=assay,
                                         pmid=pmid,
                                         SL=True)
        return sli

    def _add_reid_2016(self):
        """
         PLK1 NCBI Gene id 5347 and CKS1B NCBI Gene id 1163 
         """
        sli = SyntheticLethalInteraction(gene_A_symbol='PLK1',
                                         gene_A_id=5347,
                                         gene_B_symbol="CKS1B",
                                         gene_B_id=1163,
                                         gene_A_pert='shRNA',
                                         gene_B_pert='increased.expression',
                                         effect_type=SlConstants.GROWTH_INHIBITION_ASSAY,
                                         effect_size='n/a',
                                         cell_line='multiple.breast.cancer.cell.lines',
                                         cellosaurus_id='n/a',
                                         cancer_type=breastCarcinoma,
                                         ncit_id=ncitBreastCarcinoma,
                                         assay=SlConstants.GROWTH_INHIBITION_ASSAY,
                                         pmid='PMID:27558135',
                                         SL=True)
        self.entries.append(sli)

    def _add_parthak_2015(self):
        """
        # just 2 SL Interactions, hardcoded
        # SRC Gene is proto-oncogene blocked by Dasatinib
        # trying to maximise the Dasatinib sensitivity by SL interaction
        """
        gene1_symbol = 'SRC'
        gene1_id = 'NCBIGene:6714'
        gene1_perturbation = SlConstants.PHARMACEUTICAL
        gene2_perturbation = SlConstants.COHORT_STUDY
        pmid = 'PMID:26437225'
        assay = SlConstants.PHARAMACEUTICAL_INHIBITION_ASSAY
        effect_type = "correlation"
        cell_line = "n/a"
        cellosaurus = "n/a"
        cancer = "Recurrent Ovarian Carcinoma"
        ncit = "NCIT:C7833"
        sli = SyntheticLethalInteraction(gene_A_symbol=gene1_symbol,
                                         gene_A_id=gene1_id,
                                         gene_B_symbol="CSNK2A1",
                                         gene_B_id="NCBIGene:1457",
                                         gene_A_pert=gene1_perturbation,
                                         gene_B_pert=gene2_perturbation,
                                         effect_type=effect_type,
                                         effect_size=-0.82,
                                         cell_line=cell_line,
                                         cellosaurus_id=cellosaurus,
                                         cancer_type=cancer,
                                         ncit_id=ncit,
                                         assay=assay,
                                         pmid=pmid,
                                         SL=True)
        self.entries.append(sli)
        sli = SyntheticLethalInteraction(gene_A_symbol=gene1_symbol,
                                         gene_A_id=gene1_id,
                                         gene_B_symbol="PRKCE",
                                         gene_B_id="NCBIGene:5581",
                                         gene_A_pert=gene1_perturbation,
                                         gene_B_pert=gene2_perturbation,
                                         effect_type=effect_type,
                                         effect_size=-0.96,
                                         cell_line=cell_line,
                                         cellosaurus_id=cellosaurus,
                                         cancer_type=cancer,
                                         ncit_id=ncit,
                                         assay=assay,
                                         pmid=pmid,
                                         SL=True)
        self.entries.append(sli)


    def _add_sultana_2013(self):
        atr = 'ATR'
        atr_id = 'NCBIGene:545'
        atr_perturbation = SlConstants.PHARMACEUTICAL
        xrcc1 = 'XRCC1'
        xrcc1_id = 'NCBIGene:7515'
        xrcc1_perturbation = SlConstants.SI_RNA
        pmid = '23451157'
        assay = SlConstants.CISPLATIN_CYTOTOXICITY_ASSAY
        effect_type = "n/a"
        cell_line = "OVCAR-3"
        cellosaurus = "CVCL_0465"
        cancer = "Ovarian serous adenocarcinoma"
        ncit = "NCIT:C105555"
        sli = SyntheticLethalInteraction(gene_A_symbol=atr,
                                         gene_A_id=atr_id,
                                         gene_B_symbol=xrcc1,
                                         gene_B_id=xrcc1_id,
                                         gene_A_pert=atr_perturbation,
                                         gene_B_pert=xrcc1_perturbation,
                                         effect_type=effect_type,
                                         effect_size='n/a',
                                         cell_line=cell_line,
                                         cellosaurus_id=cellosaurus,
                                         cancer_type=cancer,
                                         ncit_id=ncit,
                                         assay=assay,
                                         pmid=pmid,
                                         SL=True)
        self.entries.append(sli)

    def _add_hoffmann_2014(self):
        """
        BRG1 and BRM identified as a SL pair in various different cancer cells
        BRM and BRG1 are closely related paralogs that function as mutually exclusive ATP-dependent catalytic subunits
        of the mSWI/SNF complexes, including A-549 (CVCL_0023)
        The current symbol for BRM is SMARCA2
        The current symbol for BRG1 is SMARCA4

        """
        smarca2 = 'SMARCA2'
        smarca2_id = 'NCBIGene:6595'
        smarca2_perturbation = SlConstants.SI_RNA
        smarca4 = 'SMARCA4'
        smarca4_id = 'NCBIGene:6597'
        smarca4_perturbation = SlConstants.LOF_MUTATION
        pmid = '31427792'
        assay = SlConstants.GROWTH_INHIBITION_ASSAY
        effect_type = "n/a"
        cell_line = "A-549"
        cellosaurus = "CVCL_0023"
        cancer = "Lung adenocarcinoma"
        ncit = "NCIT:C3512"
        sli = SyntheticLethalInteraction(gene_A_symbol=smarca2,
                                         gene_A_id=smarca2_id,
                                         gene_B_symbol=smarca4,
                                         gene_B_id=smarca4_id,
                                         gene_A_pert=smarca2_perturbation,
                                         gene_B_pert=smarca4_perturbation,
                                         effect_type=effect_type,
                                         effect_size='n/a',
                                         cell_line=cell_line,
                                         cellosaurus_id=cellosaurus,
                                         cancer_type=cancer,
                                         ncit_id=ncit,
                                         assay=assay,
                                         pmid=pmid,
                                         SL=True)
        self.entries.append(sli)


    def _add_barbie_2009(self):
        """
        Paper is about one SLI
        """
        kras = 'KRAS'
        tbk1 = 'TBK1'
        tbk1_id = 'NCBIGene:29110'
        pmid = '19847166'
        sli = self.create_sli(geneA=kras, geneAid=SlConstants.KRAS_GENE_ID, geneB=tbk1, geneBid=tbk1_id,
                              geneApert=SlConstants.ACTIVATING_MUTATION, geneBpert=SlConstants.SI_RNA,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=SlConstants.N_A, cellosaurus=SlConstants.N_A,
                              assay=SlConstants.RNA_INTERFERENCE_ASSAY,
                              cancer=SlConstants.LUNG_ADENOCARCINOMA, ncit=SlConstants.LUNG_ADENOCARCINOMA_NCIT,
                              pmid=pmid)
        self.entries.append(sli)

    def _add_scholl_2009(self):
        """
        KRAS and  STK33
        Scholl C, Fröhling S, Dunn IF, et al. Synthetic lethal interaction between oncogenic KRAS dependency and
        STK33 suppression in human cancer cells. Cell. 2009;137(5):821-834.
        """
        kras = 'KRAS'
        kras_id = 'NCBIGene:3845'
        stk33 = 'STK33'
        stk33_id = 'NCBIGene:65975'
        kras_perturbation = SlConstants.ACTIVATING_MUTATION
        stk33_perturbation = SlConstants.SI_RNA
        pmid = '19490892'
        assay = SlConstants.GROWTH_INHIBITION_ASSAY
        effect_type = SlConstants.N_A
        cell_line = SlConstants.HCT_116
        cellosaurus = SlConstants.HCT_116_CELLOSAURUS
        cancer = SlConstants.N_A
        ncit = SlConstants.N_A
        sli = SyntheticLethalInteraction(gene_A_symbol=kras,
                                         gene_A_id=kras_id,
                                         gene_B_symbol=stk33,
                                         gene_B_id=stk33_id,
                                         gene_A_pert=kras_perturbation,
                                         gene_B_pert=stk33_perturbation,
                                         effect_type=effect_type,
                                         effect_size=SlConstants.N_A,
                                         cell_line=cell_line,
                                         cellosaurus_id=cellosaurus,
                                         cancer_type=cancer,
                                         ncit_id=ncit,
                                         assay=assay,
                                         pmid=pmid,
                                         SL=True)
        self.entries.append(sli)

    def _add_chen_2020(self):
        """

        """
        pten = 'PTEN'
        pten_id = SlConstants.PTEN_GENE_ID
        mcl1 = 'MCL1'
        mcl1_id = 'NCBIGene:4170'
        pten_pertubration = SlConstants.LOF_MUTATION
        mcl1_perturbation = SlConstants.PHARMACEUTICAL
        assay = SlConstants.PHARAMACEUTICAL_INHIBITION_ASSAY
        effect_type = SlConstants.N_A
        pmid = '32737157'
        cell_line = 'isogeneic GBM cell lines'
        cellosaurus = SlConstants.N_A
        cancer = SlConstants.N_A
        ncit = SlConstants.N_A
        sli = SyntheticLethalInteraction(gene_A_symbol=pten,
                                         gene_A_id=pten_id,
                                         gene_B_symbol=mcl1,
                                         gene_B_id=mcl1_id,
                                         gene_A_pert=pten_pertubration,
                                         gene_B_pert=mcl1_perturbation,
                                         effect_type=effect_type,
                                         effect_size=SlConstants.N_A,
                                         cell_line=cell_line,
                                         cellosaurus_id=cellosaurus,
                                         cancer_type=cancer,
                                         ncit_id=ncit,
                                         assay=assay,
                                         pmid=pmid,
                                         SL=True)
        self.entries.append(sli)


    def _add_yamada_2020(self):
        arid1a = 'ARID1A'
        arid1a_id = SlConstants.ARID1A_GENE_ID
        ezh2 = 'EZH2'
        ezh2id = SlConstants.EZH2_GENE_ID
        pmid = '32506298'
        cell_line = 'isogeneic GBM cell lines'
        sli = self.create_sli(geneA=arid1a, geneAid=arid1a_id, geneB=ezh2, geneBid=ezh2id,
                              geneApert=SlConstants.LOF_MUTATION,geneBpert=SlConstants.PHARMACEUTICAL,
                              effecttype=SlConstants.N_A,effectsize=SlConstants.N_A,
                              cell=cell_line, cellosaurus=SlConstants.N_A,
                              cancer=SlConstants.GASTRIC_CARCINOMA, ncit=SlConstants.GASTRIC_CARCINOMA_NCIT,
                              assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)
        self.entries.append(sli)

    def _add_lelij_2020(self):
        stag1 = 'STAG1'
        stag1id = 'NCBIGene:10274'
        stag2 = 'STAG2'
        stag2id = 'NCBIGene:10735'
        pmid = '32467316'
        sli = self.create_sli(geneA=stag1, geneAid=stag1id, geneB=stag2, geneBid=stag2id,
                              geneApert=SlConstants.CRISPR_CAS9, geneBpert=SlConstants.LOF_MUTATION,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=SlConstants.KBM7_CELL, cellosaurus=SlConstants.KBM7_CELLOSAURUS,
                              cancer=SlConstants.CML_BCRABL_POS, ncit=SlConstants.CML_BCRABL_POS_NCIT,
                              assay=SlConstants.CRISPR_CAS9_INTERFERENCE_ASSAY, pmid=pmid)
        self.entries.append(sli)

    def _add_liu_2020(self):
        """
        Silencing p300 or CBP selectively inhibited the viability of PTEN-/- HCT116 cells, but not PTEN+/+ ones,
        indicating that the inhibition of p300 and CBP was likely to mediate the synthetic lethality phenotype induced
        by AA
        """
        pten = 'PTEN'
        ep300 = 'EP300'
        pmid = "32398948"
        sli = self.create_sli(geneA=pten, geneAid=SlConstants.PTEN_GENE_ID,
                              geneB=ep300, geneBid=SlConstants.EP300_GENE_ID,
                              geneApert=SlConstants.LOF_MUTATION, geneBpert=SlConstants.SI_RNA,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=SlConstants.HCT_116, cellosaurus=SlConstants.HCT_116_CELLOSAURUS,
                              cancer=SlConstants.COLON_CARCINOMA, ncit=SlConstants.COLON_CARCINOMA_NCIT,
                              assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)
        self.entries.append(sli)
        crebbp = 'CREBBP'
        sli = self.create_sli(geneA=crebbp, geneAid=SlConstants.CREBBP_GENE_ID,
                              geneB=ep300, geneBid=SlConstants.EP300_GENE_ID,
                              geneApert=SlConstants.LOF_MUTATION, geneBpert=SlConstants.SI_RNA,
                              effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                              cell=SlConstants.HCT_116, cellosaurus=SlConstants.HCT_116_CELLOSAURUS,
                              cancer=SlConstants.COLON_CARCINOMA, ncit=SlConstants.COLON_CARCINOMA_NCIT,
                              assay=SlConstants.CELL_VIABILITY_ASSAY, pmid=pmid)
        self.entries.append(sli)





    def get_entries(self):
        return self.entries
