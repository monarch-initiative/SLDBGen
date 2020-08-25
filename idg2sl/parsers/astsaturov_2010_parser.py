from collections import defaultdict
from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from idg2sl.gene_pair import GenePair
from .sl_constants import SlConstants
import csv


class Astsaturov2010Parser(SL_DatasetParser):
    """
    Astsaturov and colleagues identified SLs in human by combining both computational and experimental approaches.
    They first combined pathway maps, protein-protein interactions, gene expression data and human orthologs of
    Drosophila Egfr genetic interaction partners, to predict 2689 SL candidates of EGFR. They then selected 683
    candidates for RNAi screening, by their appearance in at least two of these information sources or by prior
    biological knowledge, which resulted in 61 SLs.


    """

    def __init__(self, fname=None):
        pmid = '20858866'
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        # The following genes are taken from Table 1 of the paper.
        # Validated EGFR-sensitizing genes.
        sli_list = []
        egfr = 'EGFR'
        egfr_id = self.get_ncbigene_curie(egfr)
        # Note -- I removed LOC63920 and LOC284393 from this list, I could not identify them in HGNC
        table1_symbols = {'ABL1', 'AKT2', 'ANXA6', 'ARF4', 'ARF5', 'ASCL2', 'BCAR1', 'CALM1', 'CBLC', 'CCND1', 'CD59',
                          'CDH3', 'CXCL12', 'DCN', 'DDR2', 'DIXDC1', 'DLG4', 'DUSP4', 'DUSP6', 'DUSP7', 'EPHA5',
                          'ERBB3', 'FER',  'FGFR2', 'FLNA', 'GRB7', 'HSPA9', 'INPPL1', 'KLF10',
                          'LTK', 'MAP3K1', 'MAPK1', 'MATK', 'NEDD9', 'NOTCH2', 'PIK3R1', 'PIK3R2',
                          'PIN1', 'PKN2', 'PLSCR1', 'PPIA', 'PRKACB', 'PRKCD', 'PRKCE', 'PRKCZ', 'PTPRF', 'RAC1',
                          'RAPGEF1', 'RASA3', 'RET', 'RPS6KA5', 'SC4MOL',
                          'SH2D3C', 'SHC1', 'SMAD2', 'SOS2', 'STAT3', 'TBL1Y', 'VAV3'}
        for geneB in table1_symbols:
            geneB = self.get_current_symbol(geneB)
            if geneB in self.entrez_dict:
                geneB_id = self.get_ncbigene_curie(geneB)
            else:
                raise ValueError("Could not get NCBI id for gene \"%s\" in Blomen 2015" % geneB)
            sli = SyntheticLethalInteraction(gene_A_symbol=egfr,
                                             gene_A_id=egfr_id,
                                             gene_B_symbol=geneB,
                                             gene_B_id=geneB_id,
                                             gene_A_pert=SlConstants.INHIBITORY_ANTIBODY,
                                             gene_B_pert=SlConstants.SI_RNA,
                                             effect_type=SlConstants.N_A,
                                             effect_size=0,
                                             cell_line=SlConstants.A431_CELL,
                                             cellosaurus_id=SlConstants.A431_CELLOSAURUS,
                                             cancer_type=SlConstants.N_A,
                                             ncit_id=SlConstants.N_A,
                                             assay=SlConstants.CELL_VIABILITY_ASSAY,
                                             pmid=self.pmid,
                                             SL=True)
            sli_list.append(sli)
        # The authors have one additional SLI
        # Analysis based on the Chou-Talalay coefficient of interaction showed that the small-molecule AURKA inhibitor
        # PHA-680632 (29) synergized with erlotinib in reducing cell viability of both A431 and HCT116 cells (Fig. 6B).
        # In HCT116 cells, we found strong synergy (coefficient of interaction values <0.5) between cetuximab and either
        # PHA-680632 or another AURKA inhibitor, C1368
        aurka = 'AURKA'
        aurka_id = self.get_ncbigene_curie(aurka)
        sli = SyntheticLethalInteraction(gene_A_symbol=egfr,
                                         gene_A_id=egfr_id,
                                         gene_B_symbol=aurka,
                                         gene_B_id=aurka_id,
                                         gene_A_pert=SlConstants.PHARMACEUTICAL,
                                         gene_B_pert=SlConstants.PHARMACEUTICAL,
                                         effect_type=SlConstants.N_A,
                                         effect_size=0,
                                         cell_line=SlConstants.HCT_116,
                                         cellosaurus_id=SlConstants.HCT_116_CELLOSAURUS,
                                         cancer_type=SlConstants.N_A,
                                         ncit_id=SlConstants.N_A,
                                         assay=SlConstants.CELL_VIABILITY_ASSAY,
                                         pmid=self.pmid,
                                         SL=True)
        sli_list.append(sli)
        return sli_list
