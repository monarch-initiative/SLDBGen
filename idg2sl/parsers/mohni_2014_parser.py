from collections import defaultdict
from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants
from idg2sl.gene_pair import GenePair
import csv
import numpy as np


class Mohni2014Parser(SL_DatasetParser):
    """
    The custom siRNA library targets 240 known DNA replication and DNA repair genes with four unique siRNAs per gene
    in individual wells. The authors take four experiments and say they demand that 3 of the 4 define SL. They do not
    indicate a threshold value. For our purposes, we take a threshold of -2 SD in at least 3 of 4 experiments.
    We also keep genes with a mean ABOVE zero as Negative examples.
    """

    def __init__(self, fname='data/mohniS1excerpt.tsv', entrez=None, ensembl=None, synonym=None):
        pmid = '24662920'
        super().__init__(fname=fname, pmid=pmid, entrez=entrez, ensembl=ensembl, synonym=synonym)

    def parse(self):
        geneA = 'ATR'
        geneAid = 'NCBIGene:545'
        sli_dict = defaultdict(list)
        with open(self.fname) as csvfile:
            csvreader = csv.DictReader(csvfile, delimiter='\t')
            for row in csvreader:
                geneB = row['Gene Symbol']
                geneB = self.get_current_symbol(geneB)
                if geneB == 'ATR':
                    continue # Self interaction, not a SLI!
                # A few special cases -- capitalization is not correct in the HGNC file
                if geneB == 'C10ORF119':
                    geneB = 'MCMBP'
                elif geneB == 'C15ORF20':
                    geneB = 'PIF1'
                elif geneB == 'CXORF53':
                    geneB = 'BRCC3'
                if geneB in self.entrez_dict:
                    geneBid = self.entrez_dict.get(geneB)
                else:
                    raise ValueError("Could not find id for gene %s in Mohni 2014" % geneB)
                mock1 = float(row['Mock.1'])
                atr1 = float(row['ATRi.1'])
                mock2 = float(row['Mock.2'])
                atr2 = float(row['ATRi.2'])
                mock3 = float(row['Mock.3'])
                atr3 = float(row['ATRi.3'])
                mock4 = float(row['Mock.4'])
                atr4 = float(row['ATRi.4'])
                d1 = atr1 - mock1
                d2 = atr2 - mock2
                d3 = atr3 - mock3
                d4 = atr4 - mock4
                # We demand that at least three replicates show SL
                a = np.array([d1, d2, d3, d4])
                mn = a.mean()
                if mn < -2:
                    SL = True
                elif mn >= 0:
                    SL = False
                else:
                    raise ValueError("Expecting mean either below -2 or above 0")
                sli = SyntheticLethalInteraction(gene_A_symbol=geneA,
                                                 gene_A_id=geneAid,
                                                 gene_B_symbol=geneB,
                                                 gene_B_id=geneBid,
                                                 gene_A_pert=SlConstants.PHARMACEUTICAL,
                                                 gene_B_pert=SlConstants.SI_RNA,
                                                 effect_type=SlConstants.ZSCORE,
                                                 effect_size=mn,
                                                 cell_line=SlConstants.U2OS_CELL,
                                                 cellosaurus_id=SlConstants.U2OS_CELLOSAURUS,
                                                 cancer_type='n/a',
                                                 ncit_id='n/a',
                                                 assay=SlConstants.RNA_INTERFERENCE_ASSAY,
                                                 pmid=self.pmid,
                                                 SL=SL)
                gene_pair = GenePair(geneA, geneB)
                sli_dict[gene_pair].append(sli)
        sli_list = self._mark_maximum_entries(sli_dict)
        return sli_list
