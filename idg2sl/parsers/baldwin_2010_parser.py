from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants


class Baldwin2010Parser(SL_DatasetParser):
    def __init__(self, fname=None):
        pmid = "20616055"
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        sl_genes = {'SGK2', 'PAK3'}
        tp53 = 'TP53'
        pmid = '20616055'
        tp53id = self.get_ncbigene_curie(tp53)
        sli_list = []
        for geneB in sl_genes:
            geneBid = self.get_ncbigene_curie(geneB)
            sli = SyntheticLethalInteraction(gene_A_symbol=tp53,
                                             gene_A_id=tp53id,
                                             gene_B_symbol=geneB,
                                             gene_B_id=geneBid,
                                             gene_A_pert=SlConstants.DEGRADATION,
                                             gene_B_pert=SlConstants.SH_RNA,
                                             effect_type=SlConstants.N_A,
                                             effect_size=SlConstants.N_A,
                                             cell_line='primary human foreskin keratinocytes',
                                             cellosaurus_id=SlConstants.N_A,
                                             cancer_type=SlConstants.N_A,
                                             ncit_id=SlConstants.N_A,
                                             assay=SlConstants.CELL_VIABILITY_ASSAY,
                                             pmid=pmid,
                                             SL=True)
            sli_list.append(sli)
        return sli_list