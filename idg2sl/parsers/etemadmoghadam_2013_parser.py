from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants
import csv


class Etemadmoghadam2013Parser(SL_DatasetParser):
    def __init__(self, fname=None):
        pmid = "24218601"
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        ccne1 = 'CCNE1'
        ccne1_id = self.get_ncbigene_curie(ccne1)
        sli_list = []
        sligenes = {'CDK2', 'ACAT2', 'CSE1L', 'BRCA1', 'CCNA2', 'CDC42', 'CHD2', 'DDX17', 'DUSP16', 'ENPP2',
                    'HNRNPA3', 'IARS2', 'MYC', 'PSMA5', 'RRM1', 'SLC35A3', 'SMC2', 'SPATA6', 'SRBD1', 'TPX2', 'TUBB',
                    'UBA1', 'VCP', 'XRCC2'}
        for geneB in sligenes:
            geneBId = self.get_ncbigene_curie(geneB)
            sli = SyntheticLethalInteraction(gene_A_symbol=ccne1,
                                             gene_A_id=ccne1_id,
                                             gene_B_symbol=geneB,
                                             gene_B_id=geneBId,
                                             gene_A_pert=SlConstants.OVEREXPRESSION,
                                             gene_B_pert=SlConstants.SH_RNA,
                                             effect_type=SlConstants.N_A,
                                             effect_size=SlConstants.N_A,
                                             cell_line='102 cancer cell lines',
                                             cellosaurus_id=SlConstants.N_A,
                                             cancer_type=SlConstants.N_A,
                                             ncit_id=SlConstants.N_A,
                                             assay=SlConstants.SH_RNA_DEPLETION_ASSAY,
                                             pmid=self.pmid,
                                             SL=True)
            sli_list.append(sli)
        return sli_list

