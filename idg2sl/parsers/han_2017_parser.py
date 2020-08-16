from collections import defaultdict
from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants
from idg2sl.gene_pair import GenePair
import csv


class Han2017Parser(SL_DatasetParser):
    def __init__(self, fname='data/Han2017_supplemental_table_1.tsv'):
        pmid = "28319085"
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        # using supplemental file 1
        gene1_perturbation = SlConstants.SG_RNA
        gene2_perturbation = SlConstants.SG_RNA
        assay = SlConstants.RNA_INTERFERENCE_ASSAY
        effect_type = "z-Score"
        cell_line = "K562 chronic myeloid leukemia cells"
        cellosaurus = "CVCL_0004"
        cancer = "Chronic Myelogenous Leukemia"
        ncit = "C3174"
        sli_list = []
        with open(self.fname) as csvfile:
            csvreader = csv.DictReader(csvfile, delimiter='\t')
            for row in csvreader:
                if len(row) < 4:
                    raise ValueError("Only got %d fields but was expecting at least 4 tab-separated fields" % len(row))
                # separate genes
                genes = row['Drug-target.Pairs'].split("__")
                geneA_sym = self.get_current_symbol(genes[0])
                geneB_sym = self.get_current_symbol(genes[1])
                if geneA_sym in self.entrez_dict:
                    geneA_id = "NCBIGene:{}".format(self.entrez_dict.get(geneA_sym))
                else:
                    raise ValueError("could not find id for gene A (%s) in Han 2017" % geneA_sym)
                if geneB_sym in self.entrez_dict:
                    geneB_id = "NCBIGene:{}".format(self.entrez_dict.get(geneB_sym))
                else:
                    raise ValueError("could not find id for gene B (%s) in Han 2017" % geneB_sym)
                effect = -4  # No exact value given, but authors state at least -4 for all SLIs
                sli = SyntheticLethalInteraction(gene_A_symbol=geneA_sym,
                                                 gene_A_id=geneA_id,
                                                 gene_B_symbol=geneB_sym,
                                                 gene_B_id=geneB_id,
                                                 gene_A_pert=gene1_perturbation,
                                                 gene_B_pert=gene2_perturbation,
                                                 effect_type=effect_type,
                                                 effect_size=effect,
                                                 cell_line=cell_line,
                                                 cellosaurus_id=cellosaurus,
                                                 cancer_type=cancer,
                                                 ncit_id=ncit,
                                                 assay=assay,
                                                 pmid=self.pmid,
                                                 SL=True)
                sli_list.append(sli)
        return sli_list
