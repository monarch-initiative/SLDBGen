from collections import defaultdict
from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants
from idg2sl.gene_pair import GenePair
import csv


class Kessler2012Parser(SL_DatasetParser):
    """
    We identified 403 MySL shRNAs exhibiting >2-fold decrease in abundance in the Myc-ON state
    (relative to the Myc-OFF state) (p<0.02; Fig. 1B, fig. S3, table S1).
    Note that all of the entries in the table are positives. The column 'median.pair.diffs'
    provides the LOG2 differences.
    """

    def __init__(self, fname='data/kessler2012SupplTable1.tsv'):
        pmid = 'PMID:22157079'
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        myc = 'MYC'
        myc_id = self.entrez_dict.get(myc)
        myc_perturbation = SlConstants.OVEREXPRESSION.to_string()
        geneB_perturbation = SlConstants.SI_RNA.to_string()
        assay_string = SlConstants.RNA_INTERFERENCE_ASSAY.to_string()
        effect_type = SlConstants.LOG2_DECREASE_IN_ABUNDANCE.to_string()
        cell_line = 'human mammary epithelial cells'
        cellosaurus = SlConstants.N_A.to_string()
        cancer = SlConstants.N_A.to_string()
        ncit = SlConstants.N_A.to_string()
        sli_dict = defaultdict(list)
        # Pseudogenes, divergent nc transcripts
        # DIP maps to two newer symbols (also GIF
        unclear_gene_symbols = {'ATP5EP1', 'C10orf111', 'C19ORF30', 'C3ORF51', 'CG030', 'CLEC4GP1', 'CSN1S2A', 'DIP',
                                'DKFZP434I0714', 'DVL1L1', 'FLJ20674', 'FLJ22447', 'GIF', 'HCG27', 'HMG14P',
                                'IGLV@', 'LDHBP', 'OR5D2P', 'RBMXP1', 'RPL19P1'}
        with open(self.fname) as csvfile:
            csvreader = csv.DictReader(csvfile, delimiter='\t')
            for row in csvreader:
                if len(row) != 3:
                    raise ValueError("Bad row with %d fields: %s" % (len(row), row))
                geneBsym = row['symbol']
                geneBsym = self.get_current_symbol(geneBsym)
                if geneBsym in self.entrez_dict:
                    geneB_id = "NCBIGene:{}".format(self.entrez_dict.get(geneBsym))
                elif geneBsym in unclear_gene_symbols:
                    continue
                else:
                    raise ValueError("Could not find id for %s in Kessler 2012 " % geneBsym)
                medianDiffs = float(row['median.pair.diffs'])
                sli = SyntheticLethalInteraction(gene_A_symbol=myc,
                                                 gene_A_id=myc_id,
                                                 gene_B_symbol=geneBsym,
                                                 gene_B_id=geneB_id,
                                                 gene_A_pert=myc_perturbation,
                                                 gene_B_pert=geneB_perturbation,
                                                 effect_type=effect_type,
                                                 effect_size=medianDiffs,
                                                 cell_line=cell_line,
                                                 cellosaurus_id=cellosaurus,
                                                 cancer_type=cancer,
                                                 ncit_id=ncit,
                                                 assay=assay_string,
                                                 pmid=self.pmid,
                                                 SL=True)
                gene_pair = GenePair(myc, geneBsym)
                sli_dict[gene_pair].append(sli)
        sli_list = self._mark_maximum_entries(sli_dict)
        return sli_list
