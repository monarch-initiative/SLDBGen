from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants


class Wang2019Parser(SL_DatasetParser):
    def __init__(self, fname=None):
        """

        """
        pmid = "30532030"
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        """
        overlapping identified ATRi co-essential genes. (significant in all three screens)
        RNASEH2 was validated in detail.
        Note, we replaced C17orf53 by HROB (gene id: 78995)
        MGEA5 by OGA (Gene ID: 10724)
        """

        sig_genes = {'RNASEH2B', 'RNASEH2A', "DSCC1", "TMEM208",
                     "POLE3", "POLE4", "LEO1", "CNOT1",
                     "SETD1A", "HROB", "OGA", "MCM9",
                     "USP37", "THRAP3", "DPYS", "CKS2",
                     "RHNO1", "HUS1"}
        sli_list = []
        atr = 'ATR'
        atr_id = self.get_ncbigene_curie(atr)
        for geneB in sig_genes:
            if geneB in self.entrez_dict:
                geneb_id = self.get_ncbigene_curie(geneB)
                sli = SyntheticLethalInteraction(gene_A_symbol=atr,
                                                 gene_A_id=atr_id,
                                                 gene_B_symbol=geneB,
                                                 gene_B_id=geneb_id,
                                                 gene_A_pert=SlConstants.PHARMACEUTICAL,
                                                 gene_B_pert=SlConstants.CRISPR_CAS9,
                                                 effect_type=SlConstants.N_A,
                                                 effect_size=SlConstants.N_A,
                                                 cell_line=SlConstants.N_A,
                                                 cellosaurus_id=SlConstants.N_A,
                                                 cancer_type=SlConstants.N_A,
                                                 ncit_id=SlConstants.N_A,
                                                 assay=SlConstants.CRISPR_CAS9_INTERFERENCE_ASSAY,
                                                 pmid=self.pmid,
                                                 SL=True)
                sli_list.append(sli)
            else:
                raise ValueError("Could not find id for ", geneB)
        return sli_list