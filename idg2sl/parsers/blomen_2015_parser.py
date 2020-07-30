from collections import defaultdict
from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from idg2sl.gene_pair import GenePair


class Blomen2015Parser(SL_DatasetParser):
    """
    Parse SL data from Blomen VA, , et al.
    Gene essentiality and synthetic lethality in haploid human cells.
    Science. 2015;350(6264):1092-1096.
    PMID:26472760
    synthetic lethal/sick interactions that were identified by scoring genes for fitness
    reduction in three nuclease-generated knockout clones per genotype.
    Table S7 presents synthetic lethality pairs with relatively sparse
    additional data
    Genes identified in the synthetic lethality network were analyzed for their function.
    This table lists genes that could be linked to the secretory pathway with a reference
    to the published literature or online databases.
    """
    def __init__(self, symbol2entrezID, fname, hasHeader=True):
        pmid = 'PMID:26472760'
        super().__init__(fname=fname, pmid=pmid, symbol2entrezID=symbol2entrezID, hasHeader=hasHeader)

    def parse(self):
        print("Number of rows ", len(self.rows))
        print(type(self.rows))
        perturbation = 'knockout'
        cellosuarus = 'CVCL_Y019'
        assay = 'proportions.of.sense.and.antisense.insertions'
        sli_dict = defaultdict(list)
        for row in self.rows:
            geneA = row[0]
            if geneA in self.symbol2entrezID:
                geneA_id = "NCBIGene:{}".format(self.symbol2entrezID.get(geneA))
            else:
                print("[ERROR] We could not find a gene id for " + geneA)
            geneBlist = row[3]
            for geneB in geneBlist.split(';'):
                if geneB in self.symbol2entrezID:
                    geneB_id = "NCBIGene:{}".format(self.symbol2entrezID.get(geneB))
                else:
                    print("[ERROR] We could not find a gene id for " + geneA)
                sli = SyntheticLethalInteraction(gene_A_symbol=geneA,
                                                 gene_A_id=geneA_id,
                                                 gene_B_symbol=geneB,
                                                 gene_B_id=geneB_id,
                                                 gene_A_pert=perturbation,
                                                 gene_B_pert=perturbation,
                                                 effect_type='n/a',
                                                 effect_size=0,
                                                 cell_line='HAP1',
                                                 cellosaurus_id=cellosuarus,
                                                 cancer_type='n/a',
                                                 ncit_id='n/a',
                                                 assay=assay,
                                                 pmid=self.pmid,
                                                 SL=True)
                gene_pair = GenePair(geneA, geneB)
                sli_dict[gene_pair].append(sli)
                sli_list = self._mark_maximum_entries(sli_dict)
        return sli_list
