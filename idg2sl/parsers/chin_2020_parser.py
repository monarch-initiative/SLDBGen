from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants
import csv


class Chin2020Parser(SL_DatasetParser):
    def __init__(self, fname=None):
        pmid = "33284104"
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        """
        Wnt agonist LY2090314, which mimics Wnt activation by inhibiting GSK3-β (Atkinson et al., 2015), emerged as a
        novel class of compound that inhibited the growth of all three cohesin mutants tested.
        RAD21, SMC3, and STAG2 deletion mutations in the breast epithelial cell line MCF10A resulte
        """
        gsk3b = 'GSK3B'
        gsk3b_id = self.get_ncbigene_curie(gsk3b)
        sli_list = []
        sligenes = {'RAD21', 'SMC3', 'STAG2'}
        for geneB in sligenes:
            geneBId = self.get_ncbigene_curie(geneB)
            sli = SyntheticLethalInteraction(gene_A_symbol=gsk3b,
                                             gene_A_id=gsk3b_id,
                                             gene_B_symbol=geneB,
                                             gene_B_id=geneBId,
                                             gene_A_pert=SlConstants.PHARMACEUTICAL,
                                             gene_B_pert=SlConstants.CRISPR_CAS9,
                                             effect_type=SlConstants.N_A,
                                             effect_size=SlConstants.N_A,
                                             cell_line=SlConstants.MCF10A_CELL,
                                             cellosaurus_id=SlConstants.MCF10A_CELLOSAURUS,
                                             cancer_type=SlConstants.N_A,
                                             ncit_id=SlConstants.N_A,
                                             assay=SlConstants.CELL_VIABILITY_ASSAY,
                                             pmid=self.pmid,
                                             SL=True)
            sli_list.append(sli)
        return sli_list

