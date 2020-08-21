from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants


class Krastev2011Parser(SL_DatasetParser):
    def __init__(self, fname=None):
        """
        examination of the screening results also revealed three genes that increased the wild-type/knockout ratio.
        Knockdown of UNRIP (also known asSTRAP, serine/threonine kinase receptor associated protein),
        MASTL (microtubule-associatedserine/threonine kinase-like) and KIAA1344 (also known asTXNDC16,thioredoxin
        domain containing 16) had no effect on the growth ofthe wild-type cells, whereas expression of these genes was
        required for efficient growth of the knockout cells (Fig. 2a), raising the possibility that they act in
        conjunction with TP53 to maintain cellular fitness.
        """
        pmid = "21642980"
        super().__init__(fname=fname, pmid=pmid)

    def create_sli(self, geneB_sym):
        tp53 = 'TP53'
        tp53id = self.entrez_dict.get(tp53)
        geneB_id = self.entrez_dict.get(geneB_sym)
        tp53_perturbation = SlConstants.LOF_MUTATION
        gene_B_pert = SlConstants.SI_RNA
        cell_line = SlConstants.HCT_116
        cellosaurus_id = SlConstants.HCT_116_CELLOSAURUS
        cancer_type = SlConstants.COLON_CARCINOMA
        ncit_id = SlConstants.COLON_CARCINOMA_NCIT
        assay = SlConstants.RNA_INTERFERENCE_ASSAY
        sli = SyntheticLethalInteraction(gene_A_symbol=tp53,
                                         species_id="10090",
                                         gene_A_id=tp53id,
                                         gene_B_symbol=geneB_sym,
                                         gene_B_id=geneB_id,
                                         gene_A_pert=tp53_perturbation,
                                         gene_B_pert=gene_B_pert,
                                         effect_type="n/a",
                                         effect_size="n/a",
                                         cell_line=cell_line,
                                         cellosaurus_id=cellosaurus_id,
                                         cancer_type=cancer_type,
                                         ncit_id=ncit_id,
                                         assay=assay,
                                         pmid=self.pmid,
                                         SL=True)
        return sli

    def parse(self):
        sli_list = []
        strap = 'STRAP'  # current Gene symbol for UNRIP
        sli = self.create_sli(strap)
        sli_list.append(sli)
        mastl = 'MASTL'
        sli = self.create_sli(mastl)
        sli_list.append(sli)
        txndc16 = 'TXNDC16'  # current symbol for KIAA1344
        sli = self.create_sli(txndc16)
        sli_list.append(sli)
        return sli_list
