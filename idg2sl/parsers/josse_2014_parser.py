from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants
import csv


class Josse2014Parser(SL_DatasetParser):
    """
    The authors present a screen for SL with TOP1 inhibition
    We take genes that show no evidence of this as NEGATIVEs.
    """

    def __init__(self, fname='data/josse_2014-supplement-1.tsv'):
        pmid = "25269479"
        super().__init__(fname=fname, pmid=pmid)
        self.sli_list = []

    def _add_negatives(self):
        top1 = 'TOP1'
        top1_id = self.get_ncbigene_curie(top1)
        # header Symbol	Gene_ID	Rank	RSA p-value	FDR
        with open(self.fname) as csvfile:
            csvreader = csv.DictReader(csvfile, delimiter='\t')
            for row in csvreader:
                if len(row) != 5:
                    raise ValueError("Bad line with %d instead of 6 fields: %s" % (len(row), row))
                geneB = row['Symbol']
                pval = float(row['RSA.p-value'])
                if pval < 0.5:
                    continue
                sym = self.get_current_symbol(geneB)
                if sym in self.entrez_dict:
                    # We skip symbols that cannot be identified for this negative list
                    geneB_id = self.get_ncbigene_curie(sym)
                    if top1_id == geneB_id:
                        continue  # There is one self-loop in the data, we discard it because self-loops
                        # cannot be SLIs
                    sli = SyntheticLethalInteraction(gene_A_symbol=top1,
                                                     gene_A_id=top1_id,
                                                     gene_B_symbol=sym,
                                                     gene_B_id=geneB_id,
                                                     gene_A_pert=SlConstants.PHARMACEUTICAL,
                                                     gene_B_pert=SlConstants.SI_RNA,
                                                     effect_type=SlConstants.PVAL,
                                                     effect_size=pval,
                                                     cell_line=SlConstants.MDAMB231_CELL,
                                                     cellosaurus_id=SlConstants.MDAMB231_CELLOSAURUS,
                                                     cancer_type=SlConstants.N_A,
                                                     ncit_id=SlConstants.N_A,
                                                     assay=SlConstants.CELL_VIABILITY_ASSAY,
                                                     pmid=self.pmid,
                                                     SL=False)
                    self.sli_list.append(sli)

    def parse(self):
        """
        For positives, we take genes with more than half of ≥ 7 siRNAs yielding > 4-fold sensitization (Figure 1b)
        """
        self._add_negatives()
        # TAB2 -- current symbol for MAP3K7IP2
        positive_sl = {'ATR', 'TAB2', 'PPP2R1A', 'RNF31', 'TRAF6', 'UPF1', 'USP5'}
        top1 = 'TOP1'
        top1_id = self.get_ncbigene_curie(top1)
        for geneB in positive_sl:
            geneB_id = self.get_ncbigene_curie(geneB)
            sli = SyntheticLethalInteraction(gene_A_symbol=top1,
                                             gene_A_id=top1_id,
                                             gene_B_symbol=geneB,
                                             gene_B_id=geneB_id,
                                             gene_A_pert=SlConstants.PHARMACEUTICAL,
                                             gene_B_pert=SlConstants.SI_RNA,
                                             effect_type=SlConstants.N_A,
                                             effect_size=SlConstants.N_A,
                                             cell_line=SlConstants.MDAMB231_CELL,
                                             cellosaurus_id=SlConstants.MDAMB231_CELLOSAURUS,
                                             cancer_type=SlConstants.N_A,
                                             ncit_id=SlConstants.N_A,
                                             assay=SlConstants.CELL_VIABILITY_ASSAY,
                                             pmid=self.pmid,
                                             SL=True)
            self.sli_list.append(sli)
        return self.sli_list

