from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants
import csv


class Sullivan2012Parser(SL_DatasetParser):
    """
    Nutlins are cis-imidazoline analogs which inhibit the interaction between mdm2 and tumor suppressor p53.
    Nutlin small molecules occupy p53 binding pocket of MDM2 and effectively disrupt the p53–MDM2 interaction that
    leads to activation of the p53 pathway in p53 wild-type cells. Nutlin-3 has been shown to affect the production
    of p53 within minutes. (Wikipedia).
    SLNs (read counts in DMSO > Nutlin-3).
    We found that 18/30 candidate SLNs displayed <50% relative viability and 24/30 showed <75% relative viability.
    I tried to get negative examples, but the average p value for HCT116 was p=0.366393 and strangely it was
    p=0.009691 for A549
    """

    def __init__(self, fname=None):
        pmid = "22660439"
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        sli_list = []
        tp53 = 'TP53'
        tp53id = self.get_ncbigene_curie(tp53)
        # Here we take the that 18/30 candidate SLNs displayed <50% relative viability
        # correct symbol for DICER is DICER1
        pos_sli = {'AMFR', 'ATM', 'CAPN9', 'DICER1',  'MACF1', 'MADCAM1', 'MCL1', 'MED21', 'MET',
                   'MON1B', 'PLCB4', 'RAB8B', 'RAD1', 'SRPK1', 'STAU1', 'TGFB2', 'TRPC1', 'VEGFA'}
        for geneB in pos_sli:
            genebid = self.get_ncbigene_curie(geneB)
            sli = SyntheticLethalInteraction(gene_A_symbol=tp53,
                                             gene_A_id=tp53id,
                                             gene_B_symbol=geneB,
                                             gene_B_id=genebid,
                                             gene_A_pert=SlConstants.PHARMACEUTICAL,
                                             gene_B_pert=SlConstants.SI_RNA,
                                             effect_type=SlConstants.N_A,
                                             effect_size=SlConstants.N_A,
                                             cell_line=SlConstants.HCT_116,
                                             cellosaurus_id=SlConstants.HCT_116_CELLOSAURUS,
                                             cancer_type=SlConstants.N_A,
                                             ncit_id=SlConstants.N_A,
                                             assay=SlConstants.CELL_VIABILITY_ASSAY,
                                             pmid=self.pmid,
                                             SL=True)
            sli_list.append(sli)
        return sli_list

