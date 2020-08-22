from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants


class Dai2013Parser(SL_DatasetParser):
    """
    Using a genome-wide siRNA library screening and biologic characterization, we identified that inhibition of
    thioredoxin reductase-1 (TXNRD1), one of the key antioxidant enzymes, with siRNAs or its inhibitor, auranofin,
    sensitized NSCLC cells to MK2206 treatment in vitro and in vivo. Furthermore, we found that the synthetic lethality
    interaction between the TXNRD1 and AKT pathways occurred through the KEAP1/NRF2 cellular antioxidant pathway.
    genome-wide siRNA screening was performed with the Dharmacon siRNA siGENOME library in the NSCLC cell line HCC193.
    Using a P value of 0.007542048 and at least 50% percent of growth inhibition as the cut-off, we identified 156 hits.
    Using cell proliferation assay and siRNA knockdown, we further validated the top 10 hits, including TXNRD1, CFLAR,
    ATP5J, SLC9A10, and CLN10. The results showed that knockdown of TXNRD1, CFLAR, SCLC9A10, AHCYL1, and CLN10 sensitized
    both HCC193 and H1993 cells to 1μM MK2206
    For our purposes, we will take the top ten (validated) hits. Little information is presented about the remaining data
    """

    def __init__(self, fname=None):
        pmid = "23824739"
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        # LOC400858 relates to
        # This record has been withdrawn by NCBI. This record was defined by XM_375931.4, which is not sufficient
        # evidence to define a distinct locus.
        # LOC391269 could not be identified
        # We replaced CLN10 by CTSD (current symbol of cathepsin D)
        # we replaced ATP5J by ATP5PF (current symbol of ATP synthase peripheral stalk subunit F6)
        # We replaced SLC9A10 by SLC9C1 (current symbol of solute carrier family 9 member C1)
        # We replaced EIF3S10 by EIF3A (current symbol of eukaryotic translation initiation factor 3 subunit A)
        hits = {'TXNRD1', 'CFLAR', 'SLC9C1', 'AHCYL1', 'CTSD', 'CHML', 'ATP5PF', 'EIF3A', }
        akt1 = 'AKT1'  # AKT1 is the current symbol for AKT
        akt1_id = self.get_ncbigene_curie(akt1)
        keap1 = 'KEAP1' # Back gene (wildtype) in this experiment
        keap1_id = self.get_ncbigene_curie(keap1)
        sli_List = []
        for geneB in hits:
            geneBid = self.get_ncbigene_curie(geneB)
            sli = SyntheticLethalInteraction(gene_A_symbol=akt1, gene_A_id=akt1_id,
                                             gene_B_symbol=geneB, gene_B_id=geneBid,
                                             gene_A_pert=SlConstants.PHARMACEUTICAL,
                                             gene_B_pert=SlConstants.SH_RNA,
                                             cell_line=SlConstants.HCC193_CELL,
                                             cellosaurus_id=SlConstants.HCC193_CELLOSAURUS,
                                             cancer_type=SlConstants.NON_SMALL_CELL_LUNG,
                                             ncit_id=SlConstants.NON_SMALL_CELL_LUNG_NCIT,
                                             effect_size=SlConstants.N_A, effect_type=SlConstants.N_A,
                                             assay=SlConstants.GROWTH_INHIBITION_ASSAY,
                                             background_dependency_status=SlConstants.WILDTYPE,
                                             background_dependency_gene_symbol=keap1,
                                             background_dependency_gene_id=keap1_id,
                                             SL=True,
                                             pmid=self.pmid)
            sli_List.append(sli)
        return sli_List
