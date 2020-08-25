from idg2sl import SyntheticLethalInteraction
import csv
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants



class Kim2011Parser(SL_DatasetParser):
    """
    PX-866 and NVPBEZ-235: PI3K inhibitors
    genes whose shRNA expression in untreated cells was 2 times higher than that in PX-866–treated cells were selected
    for further analysis.  To avoid off-target effects, as shown in Table 1, we selected 30 overlapping genes among
    3 cell lines identified by the screen and showed that knockdown of these genes together with PX-866 treatment could
    lead to a lethal phenotype regardless of the genetic background of the cell lines.
    Hierarchical cluster analysis of the 30 overlapping genes revealed overexpression patterns of 15 genes in specimens
    from 116 GBM patients represented in the GSE4290 dataset, as compared with 21 nontumor tissues (Fig. 2). Therefore,
    we selected the 15 genes as a subset of our target genes shown in bold type in Table 1.
    median overall survival time of all the patients with overexpression of the 15 genes was shorter
    For secondary validation, we performed a similar shRNA library screening with the U87 cell line using a dual
    PI3K/mTOR inhibitor, NVP-BEZ235. The NVP-BEZ235 screen also identified the 15 genes set similar to the PX-866
    screen, with very consistent expression levels of the inhibited targets in the comparison analysis (Table 2).
    These findings further confirm the reproducibility of the synthetic lethality screen and suggest that the 15 genes
    identified may be potential novel targets that have synergistic effects with PI3K inhibition in GBM.
    We validated one of the targets, RACK1, in this study.
    The cell viability assay showed that RACK1 knockdown cells showed enhanced sensitivity to PX-866, confirming that
    inactivation of RACK1 sensitizes glioma cells to PX-866–induced cell death, thereby showing a strong synergism of
    inactivating PI3 kinase and RACK1

    For our purposes, we will take the 15 genes in the final screen as positives. THe supplemental material for this
    paper is no longer available (link broken).

    """

    def __init__(self, fname='data/kim-2011-table2.tsv'):
        pmid = "21430111"
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        sli_list = []
        pik3ca = 'PIK3CA'
        pik3ca_id = self.get_ncbigene_curie(pik3ca)
        with open(self.fname) as csvfile:
            csvreader = csv.DictReader(csvfile, delimiter='\t')
            for row in csvreader:
                geneBsym = self.get_current_symbol(row['Gene symbol'])
                geneBid = self.get_ncbigene_curie(geneBsym)
                px866 = float(row['PX-866'])
                nvpbez235 = float(row['NVP-BEZ235'])
                mean_fc = 0.5 * (px866+nvpbez235)
                sli = SyntheticLethalInteraction(gene_A_symbol=pik3ca,
                                                 gene_A_id=pik3ca_id,
                                                 gene_B_symbol=geneBsym,
                                                 gene_B_id=geneBid,
                                                 gene_A_pert=SlConstants.PHARMACEUTICAL,
                                                 gene_B_pert=SlConstants.SI_RNA,
                                                 effect_type=SlConstants.FOLD_CHANGE,
                                                 effect_size=mean_fc,
                                                 cell_line=SlConstants.U87WT_CELL,
                                                 cellosaurus_id=SlConstants.U87WT_CELLOSAURUS,
                                                 cancer_type=SlConstants.N_A,
                                                 ncit_id=SlConstants.N_A,
                                                 assay=SlConstants.GROWTH_INHIBITION_ASSAY,
                                                 pmid=self.pmid,
                                                 SL=True)
                sli_list.append(sli)
        return sli_list