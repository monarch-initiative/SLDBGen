from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants


class Mondal2019Parser(SL_DatasetParser):
    def __init__(self, fname=None):
        pmid = "30975996"
        super().__init__(fname=fname, pmid=pmid)

    def create_sli(self, geneB, SL):
        STAG2 = 'STAG2'
        STAG2_id = self.get_ncbigene_curie(STAG2)
        geneBid = self.get_ncbigene_curie(geneB)
        sli = SyntheticLethalInteraction(gene_A_symbol=STAG2, gene_A_id=STAG2_id, gene_B_symbol=geneB, gene_B_id=geneBid,
                                         gene_A_pert=SlConstants.LOF_MUTATION, gene_B_pert=SlConstants.SH_RNA,
                                         cell_line=SlConstants.H4_CELL, cellosaurus_id=SlConstants.H4_CELLOSAURUS,
                                         cancer_type=SlConstants.N_A, ncit_id=SlConstants.N_A,
                                         effect_size=SlConstants.N_A, effect_type=SlConstants.N_A,
                                         assay=SlConstants.CELL_VIABILITY_ASSAY, SL=SL, pmid=self.pmid)
        return sli


    def parse(self):
        sli_list = []
        sl_genes = {'STAG1', 'ATR', 'BRCA1', 'RAD51', 'XRCC5', 'PRKDC'}
        for gene in sl_genes:
            sli = self.create_sli(geneB=gene, SL=True)
            sli_list.append(sli)
        not_sl_genes = {'APEX1', 'POLB', 'ERCC1', 'ERCC4'}
        for gene in not_sl_genes:
            sli = self.create_sli(geneB=gene, SL=False)
            sli_list.append(sli)
        return sli_list









        return sli_list