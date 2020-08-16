from collections import defaultdict
from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants
from idg2sl.gene_pair import GenePair
import csv


class Srivas2016Parser(SL_DatasetParser):
    """
    # using the human SL interactions (supplemental file 4 page 2)
    # https://www.cell.com/molecular-cell/fulltext/S1097-2765(16)30280-5?innerTabgraphical_S1097276516302805=#secsectitle0105
    """
    def __init__(self, fname='data/Srivas_2016.tsv'):
        pmid = '27453043'
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        gene1_perturbation = SlConstants.PHARMACEUTICAL
        gene2_perturbation = 'natural (is a TSG)'
        assay = "pharmaceutical + siRNA"
        # The following keeps track of the current largest effect size SLI for any given gene A/gene B pair
        sli_dict = defaultdict(list)
        with open(self.fname) as csvfile:
            csvreader = csv.DictReader(csvfile, delimiter='\t')
            for row in csvreader:
                if len(row) < 4:
                    raise ValueError(
                        "Only got %d fields but was expecting at least 4 tab-separated fields" % len(fields))
                # seperate col containing multiple genes
                geneA_sym = row['geneAlist'].split(",")
                geneB_sym = row['geneB']
                geneB_sym = self.get_current_symbol(geneB_sym)
                if geneB_sym in self.get_current_symbol(geneB_sym):
                    geneB_id = "NCBIGene:{}".format(self.entrez_dict.get(geneB_sym))
                else:
                    raise ValueError("Could not find id for geneB %s in Srivasa 2016" % geneB_sym)
                effect = float(row['effect'].replace(",", "."))
                for i in geneA_sym:
                    i = self.get_current_symbol(i)
                    if i in self.entrez_dict:
                        geneA_id = "NCBIGene:{}".format(self.entrez_dict.get(i))
                    else:
                        raise ValueError("Could not find id for geneA %s in Srivasa 2016" % i)
                    sli = SyntheticLethalInteraction(gene_A_symbol=i,
                                                     gene_A_id=geneA_id,
                                                     gene_B_symbol=geneB_sym,
                                                     gene_B_id=geneB_id,
                                                     gene_A_pert=gene1_perturbation,
                                                     gene_B_pert=gene2_perturbation,
                                                     effect_type=SlConstants.ZSCORE,
                                                     effect_size=effect,
                                                     cell_line=SlConstants.HELA_CELL,
                                                     cellosaurus_id=SlConstants.HELA_CELLOSAURUS,
                                                     cancer_type=SlConstants.N_A,
                                                     ncit_id=SlConstants.N_A,
                                                     assay=assay,
                                                     pmid=self.pmid,
                                                     SL=True)
                    gene_pair = GenePair(i, geneB_sym)
                    sli_dict[gene_pair].append(sli)
        sli_list = self._mark_maximum_entries(sli_dict)
        return sli_list