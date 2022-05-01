from collections import defaultdict
from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants
from idg2sl.gene_pair import GenePair
import csv


class Shen2017Parser(SL_DatasetParser):
    """
    The authors targeted all pairs of 73 cancer genes with dual-guide RNAs in three cell lines, altogether comprising 141,912
    tests of interaction. Numerous therapeutically relevant interactions were identified and these patterns replicated
     with combinatorial drugs at 75% precision.
    """

    def __init__(self, fname='data/shen2017.tsv'):
        pmid = "28319113"
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        gene1_perturbation = SlConstants.CRISPR_CAS9
        gene2_perturbation = SlConstants.CRISPR_CAS9
        assay = SlConstants.CRISPR_CAS9_INTERFERENCE_ASSAY
        effect_type = SlConstants.ZSCORE
        # The following keeps track of the current largest effect size SLI for any given gene A/gene B pair
        sli_dict = defaultdict(list)
        with open(self.fname) as csvfile:
            csvreader = csv.DictReader(csvfile, delimiter='\t')
            for row in csvreader:
                if len(row) < 7:
                    raise ValueError("Only got %d fields but was expecting at least 7 tab-separated fields" % len(row))
                geneA_sym = row['geneA']
                geneA_sym = self.get_current_symbol(geneA_sym)
                if geneA_sym in self.entrez_dict:
                    geneA_id = "NCBIGene:{}".format(self.entrez_dict.get(geneA_sym))
                else:
                    raise ValueError("Could not get gene A in Shen 2017: %s" % geneA_sym)
                geneB_sym = row['geneB']
                geneB_sym = self.get_current_symbol(geneB_sym)
                if geneB_sym in self.entrez_dict:
                    geneB_id = "NCBIGene:{}".format(self.entrez_dict.get(geneB_sym))
                else:
                    raise ValueError("Could not get gene B in Shen 2017: %s" % geneB_sym)
                if row['Interaction_type'] == "Synthetic Lethal":
                    SL = True
                else:
                    SL = False
                cell_line_list = row['Hit_Cell_Line'].split(",")
                for cell_line in cell_line_list:
                    cell_line = cell_line.strip()
                    if cell_line == "293T":
                        cellosaurus = SlConstants.CELL_293T_CELLOSAURUS
                        effect = float(row['293T_Z'].replace(",", "."))
                    elif cell_line.upper() == "HELA":
                        cell_line = SlConstants.HELA_CELL
                        cellosaurus = SlConstants.HELA_CELLOSAURUS
                        effect = float(row['HeLa_Z'].replace(",", "."))
                    elif cell_line == "A549":
                        cell_line = SlConstants.A549_CELL
                        cellosaurus = SlConstants.A549_CELLOSAURUS
                        effect = float(row['A549_Z'].replace(",", "."))
                    else:
                        raise ValueError("Could not find cell line (\"%s\") from %s" % (cell_line,row['Hit_Cell_Line']))
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
                                                     cancer_type=SlConstants.N_A,
                                                     ncit_id=SlConstants.N_A,
                                                     assay=assay,
                                                     pmid=self.pmid,
                                                     SL=SL)
                    gene_pair = GenePair(geneA_sym, geneB_sym)
                    sli_dict[gene_pair].append(sli)
        sli_list = self._mark_maximum_entries(sli_dict)
        return sli_list
