from idg2sl.parsers.sl_constants import SlConstants
from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser



class ManualEntryOne(SL_DatasetParser):
    """
       In a number of papers, only one or a handful of synthetic lethal interactions are described.
       These are very valuable. It is easiest to enter this information by hand.
       This class is basically the same as ManualEntry but split up to limit the size of the
       entries
       """

    def __init__(self, entrez, ensembl, synonym):
        super().__init__(fname=None, pmid=None, entrez=entrez, ensembl=ensembl, synonym=synonym)
        self.entries = []
        self._add_mcmanus_2009()

    def create_and_add_sli(self, geneA, geneB, geneApert, geneBpert, assay, pmid,
                           cell=SlConstants.N_A, cellosaurus=SlConstants.N_A,
                           cancer=SlConstants.N_A, ncit=SlConstants.N_A,
                           effecttype=SlConstants.N_A, effectsize=SlConstants.N_A,
                           background_dependency_status=SlConstants.N_A,
                           background_dependency_gene_symbol=SlConstants.N_A,
                           background_dependency_gene_id=SlConstants.N_A,
                           sl=True):
        geneAid = self.get_ncbigene_curie(geneA)
        geneBid = self.get_ncbigene_curie(geneB)
        sli = SyntheticLethalInteraction(gene_A_symbol=geneA,
                                         gene_A_id=geneAid,
                                         gene_B_symbol=geneB,
                                         gene_B_id=geneBid,
                                         gene_A_pert=geneApert,
                                         gene_B_pert=geneBpert,
                                         effect_type=effecttype,
                                         effect_size=effectsize,
                                         cell_line=cell,
                                         cellosaurus_id=cellosaurus,
                                         cancer_type=cancer,
                                         ncit_id=ncit,
                                         assay=assay,
                                         background_dependency_status=background_dependency_status,
                                         background_dependency_gene_symbol=background_dependency_gene_symbol,
                                         background_dependency_gene_id=background_dependency_gene_id,
                                         pmid=pmid,
                                         SL=sl)
        self.entries.append(sli)

    def _add_mcmanus_2009(self):
        pmid = '19218431'
        rad54b = 'RAD54B'  # protein kinase C, delta
        fen1 = 'FEN1'
        # These are mouse cells but other experiments were done with human cells
        # that document sufficiently the effect

        self.create_and_add_sli(geneA=rad54b, geneB=fen1, geneApert=SlConstants.ACTIVATING_MUTATION,
                                geneBpert=SlConstants.SH_RNA, cell=SlConstants.HCT_116, cellosaurus=SlConstants.HCT_116_CELLOSAURUS,
                                assay=SlConstants.GROWTH_INHIBITION_ASSAY, pmid=pmid)


    def get_entries(self):
        return self.entries
