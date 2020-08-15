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
        self._get_reid_2016()
        self._get_parthak_2015()
        self._get_sultana_2013()

    def _get_reid_2016(self):
        """
         PLK1 NCBI Gene id 5347 and CKS1B NCBI Gene id 1163 
         """
        sli = SyntheticLethalInteraction(gene_A_symbol='PLK1',
                                         gene_A_id=5347,
                                         gene_B_symbol="CKS1B",
                                         gene_B_id=1163,
                                         gene_A_pert='shRNA',
                                         gene_B_pert='increased.expression',
                                         effect_type=SlConstants.GROWTH_INHIBITION_ASSAY.to_string(),
                                         effect_size='n/a',
                                         cell_line='multiple.breast.cancer.cell.lines',
                                         cellosaurus_id='n/a',
                                         cancer_type=breastCarcinoma,
                                         ncit_id=ncitBreastCarcinoma,
                                         assay=SlConstants.GROWTH_INHIBITION_ASSAY.to_string(),
                                         pmid='PMID:27558135',
                                         SL=True)
        self.entries.append(sli)

    def _get_parthak_2015(self):
        """
        # just 2 SL Interactions, hardcoded
        # SRC Gene is proto-oncogene blocked by Dasatinib
        # trying to maximise the Dasatinib sensitivity by SL interaction
        """
        gene1_symbol = 'SRC'
        gene1_id = 'NCBIGene:6714'
        gene1_perturbation = SlConstants.PHARMACEUTICAL.to_string()
        gene2_perturbation = SlConstants.COHORT_STUDY.to_string()
        pmid = 'PMID:26437225'
        assay = SlConstants.PHARAMACEUTICAL_INHIBITION_ASSAY.to_string()
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


    def _get_sultana_2013(self):
        atr = 'ATR'
        atr_id = 'NCBIGene:545'
        atr_perturbation = SlConstants.PHARMACEUTICAL.to_string()
        xrcc1 = 'XRCC1'
        xrcc1_id = 'NCBIGene:7515'
        xrcc1_perturbation = SlConstants.SI_RNA.to_string()
        pmid = '23451157'
        assay = SlConstants.CISPLATIN_CYTOTOXICITY_ASSAY.to_string()
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

    def _get_hoffmann_2014(self):
        """
        BRG1 and BRM identified as a SL pair in various different cancer cells
        BRM and BRG1 are closely related paralogs that function as mutually exclusive ATP-dependent catalytic subunits
        of the mSWI/SNF complexes, including A-549 (CVCL_0023)
        The current symbol for BRM is SMARCA2
        The current symbol for BRG1 is SMARCA4

        """
        smarca2 = 'SMARCA2'
        smarca2_id = 'NCBIGene:6595'
        smarca2_perturbation = SlConstants.SI_RNA.to_string()
        smarca4 = 'SMARCA4'
        smarca4_id = 'NCBIGene:6597'
        smarca4_perturbation = SlConstants.LOF_MUTATION.to_string()
        pmid = '31427792'
        assay = SlConstants.GROWTH_INHIBITION_ASSAY.to_string()
        effect_type = "n/a"
        cell_line = "A-549"
        cellosaurus = "CVCL_0023"
        cancer = "Lung adenocarcinoma"
        ncit = "NCIT:C3512"
        sli = SyntheticLethalInteraction(gene_A_symbol=smarca2,
                                         gene_A_id=smarca2_id,
                                         gene_B_symbol=smarca4,
                                         gene_B_id=smarca4,
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


    def _get_barbie_2009(self):
        """
        Paper is about one SLI
        """
        kras = 'KRAS'
        kras_id = 'NCBIGene:3845'
        tbk1 = 'TBK1'
        tbk1_id = 'NCBIGene:29110'
        kras_perturbation = SlConstants.ACTIVATING_MUTATION.to_string()
        tbk1_perturbation = SlConstants.SI_RNA.to_string()
        pmid = '19847166'
        assay = SlConstants.GROWTH_INHIBITION_ASSAY.to_string()
        effect_type = "n/a"
        cell_line = "n/a"
        cellosaurus = "n/a"
        cancer = "Lung adenocarcinoma"
        ncit = "NCIT:C3512"
        sli = SyntheticLethalInteraction(gene_A_symbol=kras,
                                         gene_A_id=kras_id,
                                         gene_B_symbol=tbk1,
                                         gene_B_id=tbk1_id,
                                         gene_A_pert=kras_perturbation,
                                         gene_B_pert=tbk1_perturbation,
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


    def get_entries(self):
        return self.entries
