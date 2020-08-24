from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants



class Wang2016Parser(SL_DatasetParser):
    """
    This paper presents mainly a computational analysis but 3 SLIs are validated convincingly.
    """

    def __init__(self, fname=None):
        pmid = "26937901"
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        sli_list = []
        apc = 'APC'
        apc_id = self.get_ncbigene_curie(apc)
        sli_genes = {'TDO2', 'CTNNB1', 'CSNK1A1'}
        for geneB in sli_genes:
            genebid = self.get_ncbigene_curie(geneB)
            sli = SyntheticLethalInteraction(gene_A_symbol=apc,
                                             gene_A_id=apc_id,
                                             gene_B_symbol=geneB,
                                             gene_B_id=genebid,
                                             gene_A_pert=SlConstants.PHARMACEUTICAL,
                                             gene_B_pert=SlConstants.SI_RNA,
                                             effect_type=SlConstants.N_A,
                                             effect_size=SlConstants.N_A,
                                             cell_line=SlConstants.SW620_CELL,
                                             cellosaurus_id=SlConstants.SW620_CELLOSAURUS,
                                             cancer_type=SlConstants.N_A,
                                             ncit_id=SlConstants.N_A,
                                             assay=SlConstants.GROWTH_INHIBITION_ASSAY,
                                             pmid=self.pmid,
                                             SL=True)
            sli_list.append(sli)
        return sli_list

