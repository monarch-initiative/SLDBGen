from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants


class Mengwasser2019Parser(SL_DatasetParser):
    def __init__(self, fname=None):
        pmid = " 30686591"
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        """
        While POLQ serves as a positive control, FEN1 and APEX2 represent novel B2SL genes and novel potential drug
        targets in BRCA-deficient tumors. The authors do not concretely name the entire list of SLIs, so
        we restrict ourselves to the three that are investigated in detail. FEN1 was also validated for BRCA1
        """
        brca2 = 'BRCA2'
        brca2_id = self.get_ncbigene_curie(brca2)
        geneBlist = {'POLQ', 'FEN1', 'APEX2'}
        sli_list = []
        for geneB in geneBlist:
            geneB_id = self.get_ncbigene_curie(geneB)
            sli = SyntheticLethalInteraction(gene_A_symbol=brca2,
                                             gene_A_id=brca2_id,
                                             gene_B_symbol=geneB,
                                             gene_B_id=geneB_id,
                                             gene_A_pert=SlConstants.LOF_MUTATION,
                                             gene_B_pert=SlConstants.CRISPR_CAS9,
                                             effect_type=SlConstants.N_A,
                                             effect_size=SlConstants.N_A,
                                             cell_line=SlConstants.PEO1_CELL,
                                             cellosaurus_id=SlConstants.PEO1_CELLOSAURUS,
                                             cancer_type=SlConstants.N_A,
                                             ncit_id=SlConstants.N_A,
                                             assay=SlConstants.MULTICOLOR_COMPETITION_ASSAY,
                                             pmid=self.pmid,
                                             SL=True)
            sli_list.append(sli)
        brca1 = 'BRCA1'
        brca1_id = self.get_ncbigene_curie(brca1)
        geneBlist = {'FEN1', 'APEX2'}
        for geneB in geneBlist:
            geneB_id = self.get_ncbigene_curie(geneB)
            sli = SyntheticLethalInteraction(gene_A_symbol=brca1,
                                         gene_A_id=brca1_id,
                                         gene_B_symbol=geneB,
                                         gene_B_id=geneB_id,
                                         gene_A_pert=SlConstants.LOF_MUTATION,
                                         gene_B_pert=SlConstants.CRISPR_CAS9,
                                         effect_type=SlConstants.N_A,
                                         effect_size=SlConstants.N_A,
                                         cell_line='BRCA1 isogenic RPE1 cell line',
                                         cellosaurus_id=SlConstants.N_A,
                                         cancer_type=SlConstants.N_A,
                                         ncit_id=SlConstants.N_A,
                                         assay=SlConstants.MULTICOLOR_COMPETITION_ASSAY,
                                         pmid=self.pmid,
                                         SL=True)
            sli_list.append(sli)
        return sli_list
