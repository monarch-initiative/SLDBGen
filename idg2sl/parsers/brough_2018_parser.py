from collections import defaultdict
from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants
from idg2sl.gene_pair import GenePair
import csv


class Brough2018Parser(SL_DatasetParser):
    """
    Using the siMEM approach with our Rb classification of 42 TNBC TCLs, we identified 1065 Rb-specific dependencies
    (p < 0.05, siMEM, Fig. 2c): 437 genes, where shRNA preferentially inhibited Rb-defective TNBC (i.e. Rb synthetic
    lethal effects) and 628 genes, where shRNAs preferentially targeted Rb-proficient TNBC TCLs (Supplementary Data 5).
    We also carried out similar analyses in siRNA/shRNA data sets that included TNBC TCLs, from other sources: the
    DRIVE data set [29]; the Achilles data set [23]; and the ICR-Intercell data set [20, 21], and provide the lists
    of Rb dependencies identified in these data sets in Supplementary Data 6, 7, 8, respectively.
    To eliminate such effects from further study, we applied a pragmatic approach that removed from further assessment
    p < 0.05 synthetic lethal effects where the median zGARP score in the Rb-defective TCLs was >−1 (i.e., effects,
    where profound cell inhibition in Rb TCLs not observed) and those effects where median zGARP score in Rb-proficient
    effects was <−2 (i.e., dependencies that still elicited profound cell inhibition in Rb-proficient TCLs);
    three examples that fulfilled these criteria, GPS1, SNRPF and SNW1, are shown in Fig. 2f, g, h. This triage step
    identified 122 Rb synthetic lethal effects in the Colt2 data set that fulfilled these criteria
     (Supplementary Data 9). Similarly, triaged dependencies were identified in the DRIVE [29], Achilles [23] and
     ICR-Intercell data sets [20, 21] (Supplementary Data 10, 11 and 12, respectively).
     We will concentrate on the highly penetrant SLIs: (>80% penetrant) from all datasets

    """

    def __init__(self, fname=None):
        pmid = 'PMID:29915391'
        super().__init__(fname=fname, pmid=pmid)
        self.sli_dict = defaultdict(list)
        # old symbols that are assigned to multiple genes
        # Pseudogenes: PI4KAP2, BTF3P10, GLRA4, RPL21P44
        self.unclear_gene_symbols = {'SARS', 'PI4KAP2', 'BTF3P10', 'DLGAP1-AS1', 'GLRA4', 'URGCP-MRPS24', 'PCDHA3',
                                     'RPL21P44', 'RPL21P28', 'ZNF733P'}

    def parse_suppl9(self):
        fname = 'data/brough_2012_suppl9.tsv'
        rb1 = 'RB1'
        rb1_id = self.entrez_dict.get(rb1)
        rb1_perturbation = SlConstants.LOF_MUTATION
        gene2_perturbation = SlConstants.SI_RNA
        assay_string = "siMEM+penetrance"
        effect_type = "penetrance"
        cell_line = SlConstants.N_A
        cellosaurus = SlConstants.N_A
        cancer = SlConstants.N_A
        ncit = SlConstants.N_A

        with open(fname) as csvfile:
            csvreader = csv.DictReader(csvfile, delimiter='\t')
            for row in csvreader:
                # print(row)
                geneBsym = self.get_current_symbol(row['symbol'])
                if geneBsym in self.entrez_dict:
                    geneB_id = "NCBIGene:{}".format(self.entrez_dict.get(geneBsym))
                elif geneBsym in self.unclear_gene_symbols:
                    continue
                else:
                    raise ValueError("Could not find iid for %s in Brough 2018 2008 " % geneBsym)
                penetrance = int(row['Penetrance.%'])
                if penetrance >= 80:
                    sli = SyntheticLethalInteraction(gene_A_symbol=rb1,
                                                     gene_A_id=rb1_id,
                                                     gene_B_symbol=geneBsym,
                                                     gene_B_id=geneB_id,
                                                     gene_A_pert=rb1_perturbation,
                                                     gene_B_pert=gene2_perturbation,
                                                     effect_type=effect_type,
                                                     effect_size=penetrance,
                                                     cell_line=cell_line,
                                                     cellosaurus_id=cellosaurus,
                                                     cancer_type=cancer,
                                                     ncit_id=ncit,
                                                     assay=assay_string,
                                                     pmid=self.pmid,
                                                     SL=True)
                    gene_pair = GenePair(rb1, geneBsym)
                    self.sli_dict[gene_pair].append(sli)

    def parse_suppl10_11(self, fname):
        rb1 = 'RB1'
        rb1_id = self.entrez_dict.get(rb1)
        rb1_perturbation = SlConstants.LOF_MUTATION
        gene2_perturbation = SlConstants.SI_RNA
        assay_string = "siMEM+penetrance"
        effect_type = "penetrance"
        cell_line = SlConstants.N_A
        cellosaurus = SlConstants.N_A
        cancer = SlConstants.N_A
        ncit = SlConstants.N_A

        with open(fname) as csvfile:
            csvreader = csv.DictReader(csvfile, delimiter='\t')
            for row in csvreader:
                # print(row)
                geneBsym = self.get_current_symbol(row['target'])
                if ',' in geneBsym:
                    continue # We cannot assign an effect unambiguously to one of the genes
                    # some of the entries are like  PMS2,PMS2CL
                if geneBsym in self.entrez_dict:
                    geneB_id = "NCBIGene:{}".format(self.entrez_dict.get(geneBsym))
                elif geneBsym in self.unclear_gene_symbols:
                    continue
                else:
                    raise ValueError("Could not find id for %s in Brough 2018 2008 " % geneBsym)
                penetrance = int(row['Penetrance.(%)'])
                if penetrance >= 80:
                    sli = SyntheticLethalInteraction(gene_A_symbol=rb1,
                                                     gene_A_id=rb1_id,
                                                     gene_B_symbol=geneBsym,
                                                     gene_B_id=geneB_id,
                                                     gene_A_pert=rb1_perturbation,
                                                     gene_B_pert=gene2_perturbation,
                                                     effect_type=effect_type,
                                                     effect_size=penetrance,
                                                     cell_line=cell_line,
                                                     cellosaurus_id=cellosaurus,
                                                     cancer_type=cancer,
                                                     ncit_id=ncit,
                                                     assay=assay_string,
                                                     pmid=self.pmid,
                                                     SL=True)
                    gene_pair = GenePair(rb1, geneBsym)
                    self.sli_dict[gene_pair].append(sli)

    def parse(self):
        self.parse_suppl9()
        self.parse_suppl10_11(fname='data/brough_2012_suppl10.tsv')
        self.parse_suppl10_11(fname='data/brough_2012_suppl11.tsv')
        sli_list = self._mark_maximum_entries(self.sli_dict)
        return sli_list
