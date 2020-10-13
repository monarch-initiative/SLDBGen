from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants


class Kang2015Parser(SL_DatasetParser):
    """
    We performed a systematic RNAi screen in two BRAF V600D/E mutant expressing melanoma cells lines including WM2664
    (V600D) and A375 (V600E), as well as control BRAF WT-expressing melanoma cells including PMWK, CHL-1, HMCB
    (NRas Q61K) and SK-MEL-2 (NRas Q61R) (Figure 1B). We then used 2 RIGER methods (RIGER_SB and RIGER_KS) and
    overlapped these against another method, Gene Set Analysis R package to analyze the normalized B-scores for each
    cell line (Barbie et al., 2009; Gould et al., 2006; Malo et al., 2006; Sims et al., 2011). The top-ranked
    100 genes identified by each method were overlapped and 36 genes were enriched as top candidate synthetic lethal
    partners of BRAF V600E (Figure 1B and Tables S3-S6). In a secondary screen, we validated the 36 candidates
     (186 shRNAs) using additional BRAF V600E expressing SK-MEL-5 and A2058 melanoma cells, compared to control BRAF
     WT-expressing PMWK and HMCB cells (Figure 1B). Results analyzed by RIGER_SB, RIGER_KS and Gene Set Analysis R
     package, and 8 genes were enriched using the top 15 genes identified in the primary screen
     We extract these 8 genes
     The authors provide additional evidence for HMGCS1 in PMID: 28468827
    """
    def __init__(self, fname=None):
        pmid = "26145173"
        super().__init__(fname=fname, pmid=pmid)

    def get_sli(self, geneB_sym, pmid):
        # Gene A should be eighter NRAS or KRAS.
        # These genes had activating mutations in the cell lines
        braf = 'BRAF'
        brafID = SlConstants.BRAF_GENE_ID
        geneB_id = self.get_ncbigene_curie(geneB_sym)
        sli = SyntheticLethalInteraction(gene_A_symbol=braf,
                                         species_id="10090",
                                         gene_A_id=brafID,
                                         gene_B_symbol=geneB_sym,
                                         gene_B_id=geneB_id,
                                         gene_A_pert=SlConstants.ACTIVATING_MUTATION,
                                         gene_B_pert=SlConstants.SH_RNA,
                                         effect_type=SlConstants.N_A,
                                         effect_size=SlConstants.N_A,
                                         cell_line=SlConstants.A375_CELL,
                                         cellosaurus_id=SlConstants.A375_CELLOSAURUS,
                                         cancer_type=SlConstants.MELANOMA,
                                         ncit_id=SlConstants.MELANOMA_NCIT,
                                         assay=SlConstants.CELL_VIABILITY_ASSAY,
                                         pmid=pmid,
                                         SL=True)
        return sli


    def parseHMGCS1(self):
        pmid = '28468827'
        geneb = 'HMGCS1'
        return self.get_sli(geneB_sym=geneb, pmid=pmid)

    def parse(self):
        sli_list = []
        sli = self.parseHMGCS1()
        sli_list.append(sli)
        pmid = '26145173'
        # genes taken from Figure 1 of PMID: 26145173
        # correct symbol for CSGLCA-T is CHPF2
        geneBset = {'HMGCL', 'CYP39A1', 'CYP2C9', 'CYP2E1', 'CYP2J2', 'CYP2S1', 'CHPF2'}
        for geneB in geneBset:
            sli = self.get_sli(geneB_sym=geneB, pmid=pmid)
            sli_list.append(sli)
        return sli_list
