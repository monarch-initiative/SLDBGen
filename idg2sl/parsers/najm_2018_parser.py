from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants
from csv import DictReader
from collections import defaultdict


class Najm2018Parser(SL_DatasetParser):
    """
    Najm FJ Orthologous CRISPR-Cas9 enzymes for combinatorial genetic screens.
    Nat Biotechnol. 2018 Feb;36(2):179-189. doi: 10.1038/nbt.4048. Epub 2017 Dec 18. PMID: 29251726; PMCID: PMC5800952.

    """

    def __init__(self, fname='data/SupplementaryTable4NajmSynLetFDRs.txt'):
        pmid = "29251726"
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        """
        While POLQ serves as a positive control, FEN1 and APEX2 represent novel B2SL genes and novel potential drug
        targets in BRCA-deficient tumors. The authors do not concretely name the entire list of SLIs, so
        we restrict ourselves to the three that are investigated in detail. FEN1 was also validated for BRCA1
        """

        sli_list = []
        seen_interactions = defaultdict(float)  # There are duplicates. Keep the one with the lowest q value
        with open(self.fname) as f:
            r = DictReader(f, delimiter='\t')
            for row in r:
                gene1 = row['Gene 1']
                gene2 = row['Gene 2']
                synletQ = float(row['SynLet q-value'])
                if synletQ > 0.05:
                    continue
                if gene1 == gene2:
                    continue
                genes = sorted([gene1, gene2])
                gene_key = f"{genes[0]}-{genes[1]}"
                if gene_key not in seen_interactions:
                    seen_interactions[gene_key] = synletQ
                else:
                    seen_interactions[gene_key] = min(synletQ, seen_interactions.get(gene_key))
                    continue
        for gene_key, qval in seen_interactions.items():
            genes = gene_key.split("-")
            geneA = genes[0]
            geneA_id = self.get_ncbigene_curie(geneA)
            geneB = genes[1]
            geneB_id = self.get_ncbigene_curie(geneB)
            sli = SyntheticLethalInteraction(gene_A_symbol=geneA,
                                             gene_A_id=geneA_id,
                                             gene_B_symbol=geneB,
                                             gene_B_id=geneB_id,
                                             gene_A_pert=SlConstants.CRISPR_CAS9,
                                             gene_B_pert=SlConstants.CRISPR_CAS9,
                                             effect_type=SlConstants.QVAL,
                                             effect_size=str(qval),
                                             cell_line=SlConstants.N_A,
                                             cellosaurus_id=SlConstants.N_A,
                                             cancer_type=SlConstants.N_A,
                                             ncit_id=SlConstants.N_A,
                                             assay=SlConstants.BIG_PAPI,
                                             pmid=self.pmid,
                                             SL=True)
            sli_list.append(sli)
        return sli_list
