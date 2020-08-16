from collections import defaultdict
from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants
from idg2sl.gene_pair import GenePair
import csv


class Shen2015Parser(SL_DatasetParser):
    """
    # GeneA is always CHEK1 (pharmaceutically inhibited by AZD7762)
    # GeneB is in data/Shen_2015.tsv
    """

    def __init__(self, fname='data/Shen_2015.tsv'):
        pmid = '26437225'
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        geneA_symbol = 'CHEK1'
        geneA_id = 'NCBIGene:1111'
        geneA_perturbation = SlConstants.PHARMACEUTICAL
        gene2_perturbation = SlConstants.SI_RNA
        assay = SlConstants.RNA_INTERFERENCE_ASSAY
        effect_type = SlConstants.ZSCORE
        cell_line = "HeLa-Cells"
        cellosaurus = "CVCL_0030"
        cancer = ""
        ncit = ""  #
        sli_list = []
        # The following keeps track of the current largest effect size SLI for any given gene A/gene B pair
        sli_dict = defaultdict(list)
        with open(self.fname) as csvfile:
            csvreader = csv.DictReader(csvfile, delimiter='\t')
            # Z-Score	Symbol	Entrez ID	Gene Name
            for row in csvreader:
                if len(row) < 3:
                    raise ValueError("Only got %d fields but was expecting at least 3" % len(row))
                geneB_sym = row['Symbol']
                geneB_sym = self.get_current_symbol(geneB_sym)
                if geneB_sym in self.entrez_dict:
                    geneB_id = "NCBIGene:{}".format(self.entrez_dict.get(geneB_sym))
                else:
                    raise ValueError("Could not find id for gene symbol %s in Shen 2015" % geneB_sym)
                effect = float(row['Z-Score'].replace(",", "."))
                sl_genes = ["FZR1", "RAD17", "RFC1", "BLM", "CDC73", "CDC6", "WEE1"]
                if geneB_sym in sl_genes:
                    SL = True
                else:
                    SL = False
                sli = SyntheticLethalInteraction(gene_A_symbol=geneA_symbol,
                                                 gene_A_id=geneA_id,
                                                 gene_B_symbol=geneB_sym,
                                                 gene_B_id=geneB_id,
                                                 gene_A_pert=geneA_perturbation,
                                                 gene_B_pert=gene2_perturbation,
                                                 effect_type=effect_type,
                                                 effect_size=effect,
                                                 cell_line=cell_line,
                                                 cellosaurus_id=cellosaurus,
                                                 cancer_type=cancer,
                                                 ncit_id=ncit,
                                                 assay=assay,
                                                 pmid=self.pmid,
                                                 SL=SL)
                gene_pair = GenePair(geneA_symbol, geneB_sym)
                sli_dict[gene_pair].append(sli)
        sli_list = self._mark_maximum_entries(sli_dict)
        return sli_list
