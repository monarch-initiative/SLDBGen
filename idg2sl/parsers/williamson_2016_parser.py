from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants
import csv


class Williamson2016Parser(SL_DatasetParser):
    """

    """
    def __init__(self, fname='data/williamson-2016-suppl2.csv'):
        pmid = "27958275"
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        """
        symbol	MCF12A.Z-score	HCC1143.Z-score
        """
        sli_list = []
        atr = 'ATR'
        atr_id = self.get_ncbigene_curie(atr)
        with open(self.fname) as csvfile:
            csvreader = csv.DictReader(csvfile, delimiter='\t')
            for row in csvreader:
                geneBsym = self.get_current_symbol(row['symbol'])
                if geneBsym == 'C9ORF96':
                    geneBsym = 'STKLD1'
                geneB_id = self.get_ncbigene_curie(geneBsym)
                mcf12 = float(row['MCF12A.Z-score'])
                hcc1143 = float(row['HCC1143.Z-score'])
                meanz = 0.5 * (mcf12 + hcc1143)
                sli = SyntheticLethalInteraction(gene_A_symbol=atr,
                                                 gene_A_id=atr_id,
                                                 gene_B_symbol=geneBsym,
                                                 gene_B_id=geneB_id,
                                                 gene_A_pert=SlConstants.PHARMACEUTICAL,
                                                 gene_B_pert=SlConstants.SI_RNA,
                                                 effect_type=SlConstants.ZSCORE,
                                                 effect_size=meanz,
                                                 cell_line=SlConstants.HCC1143_CELL,
                                                 cellosaurus_id=SlConstants.HCC1143_CELLOSAURUS,
                                                 cancer_type=SlConstants.N_A,
                                                 ncit_id=SlConstants.N_A,
                                                 assay=SlConstants.CELL_VIABILITY_ASSAY,
                                                 pmid=self.pmid,
                                                 SL=True)
                sli_list.append(sli)
            return sli_list
