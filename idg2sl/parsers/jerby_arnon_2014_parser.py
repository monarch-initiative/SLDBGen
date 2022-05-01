from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants
import csv


class JerbyArnon2014Parser(SL_DatasetParser):
    """
    Jerby-Arnon L, et al., Predicting cancer-specific vulnerability via data-driven detection of synthetic lethality.
    Cell. 2014 Aug 28;158(5):1199-1209. doi: 10.1016/j.cell.2014.07.027. PMID: 25171417.
    We use this threshold
    Selective genes are those with differential inhibition scores R4.8 (the score of the positive control
    gene MYT1, identified in Bommi-Reddy et al., 2008).
    """
    def __init__(self, fname='data/JerbyArnon2014-suppl2.csv'):
        pmid = "25171417"
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        # using supplemental file 1
        gene1_perturbation = SlConstants.LOF_MUTATION
        gene2_perturbation = SlConstants.SI_RNA
        geneA_sym = 'VHL'
        assay = SlConstants.RNA_INTERFERENCE_ASSAY
        sli_list = []
        PERCENT_INHIBTION_THRESHOLD = 4.8
        current_symbols = {'SETD8': 'KMT5A', 'ADRBK1': 'GRK2', 'BZRAP1': 'TSPOAP1',
                           'MKL1': 'MRTFA'}
        with open(self.fname) as csvfile:
            csvreader = csv.DictReader(csvfile)
            for row in csvreader:
                if len(row) < 4:
                    raise ValueError("Only got %d fields but was expecting at least 4 tab-separated fields" % len(row))
                geneB_sym = row['Current Gene Symbol']
                if geneB_sym in current_symbols:
                    geneB_sym = current_symbols.get(geneB_sym)
                if geneA_sym in self.entrez_dict:
                    geneA_id = "NCBIGene:{}".format(self.entrez_dict.get(geneA_sym))
                else:
                    raise ValueError("could not find id for gene A (%s) in Han 2017" % geneA_sym)
                if geneB_sym in self.entrez_dict:
                    geneB_id = "NCBIGene:{}".format(self.entrez_dict.get(geneB_sym))
                else:
                    raise ValueError("could not find id for gene B (%s) in Han 2017" % geneB_sym)
                percent_inhibition = float(row['percent_inhib_HIGH_CONT_Median'])
                if percent_inhibition < PERCENT_INHIBTION_THRESHOLD:
                    continue
                sli = SyntheticLethalInteraction(gene_A_symbol=geneA_sym,
                                                 gene_A_id=geneA_id,
                                                 gene_B_symbol=geneB_sym,
                                                 gene_B_id=geneB_id,
                                                 gene_A_pert=gene1_perturbation,
                                                 gene_B_pert=gene2_perturbation,
                                                 effect_type=SlConstants.PERCENT_INHIBITION,
                                                 effect_size=percent_inhibition,
                                                 cell_line=SlConstants.RCC4_CELL,
                                                 cellosaurus_id=SlConstants.RCC4_CELLOSAURUS,
                                                 cancer_type=SlConstants.RENAL_CELL_CARCINOMA,
                                                 ncit_id=SlConstants.RENAL_CELL_CARCINOMA_NCIT,
                                                 assay=assay,
                                                 pmid=self.pmid,
                                                 SL=True)
                sli_list.append(sli)
            return sli_list
