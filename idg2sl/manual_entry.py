from .parsers.sl_constants import SlConstants
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

        sli_list = []

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

    def get_entries(self):
        return self.entries
