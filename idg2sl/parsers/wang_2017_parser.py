from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants



class Wang2017Parser(SL_DatasetParser):
    def __init__(self, fname=None):
        """
        Somewhat difficult paper to analyse. The results are not always clearly presented.
        revealed five genes that were required only in the context of oncogenic Ras.
        Two genes (RCE1 and ICMT) are involved in the maturation of Ras. Two additional genes
        (RAF1 and SHOC2) are involved in MAPK pathway signaling. The final gene, PREX1, did not
        immediately fit in either category and is discussed later in its own section.
        """
        pmid="28162770"
        super().__init__(fname=fname, pmid=pmid)


    def get_sli(self, geneA_sym, geneA_id, geneB_sym, geneB_id):
        # Gene A should be eighter NRAS or KRAS.
        # These genes had activating mutations in the cell lines
        ncit = "NCIT:C3171"
        cancer = "Acute Myeloid Leukemia"
        sli = SyntheticLethalInteraction(gene_A_symbol=geneA_sym,
                                         species_id="10090",
                                         gene_A_id=geneA_id,
                                         gene_B_symbol=geneB_sym,
                                         gene_B_id=geneB_id,
                                         gene_A_pert=SlConstants.ACTIVATING_MUTATION,
                                         gene_B_pert=SlConstants.SG_RNA,
                                         effect_type="n/a",
                                         effect_size="n/a",
                                         cell_line="n/a",
                                         cellosaurus_id="n/a",
                                         cancer_type=cancer,
                                         ncit_id=ncit,
                                         assay=SlConstants.CRISPR_CAS9_INTERFERENCE_ASSAY,
                                         pmid=self.pmid,
                                         SL=True)
        return sli

    def parse(self):
        nras = "NRAS"
        nrasid = "NCBIGene:4893"
        kras = "KRAS"
        krasid = "NCBIGene:3845"
        # The following keeps track of the current largest effect size SLI for any given gene A/gene B pair
        sli_list = []

        rce1 = "RCE1"
        rce1id = "NCBIGene:{}".format(self.entrez_dict.get(rce1))

        sli = self.get_sli(geneA_sym=nras, geneA_id=nrasid, geneB_sym=rce1, geneB_id=rce1id)
        sli_list.append(sli)
        sli = self.get_sli(geneA_sym=kras, geneA_id=krasid, geneB_sym=rce1, geneB_id=rce1id)
        sli_list.append(sli)

        icmt = 'ICMT'
        icmtid = "NCBIGene:{}".format(self.entrez_dict.get(icmt))
        sli = self.get_sli(geneA_sym=nras, geneA_id=nrasid, geneB_sym=icmt, geneB_id=icmtid)
        sli_list.append(sli)
        sli = self.get_sli(geneA_sym=kras, geneA_id=krasid, geneB_sym=icmt, geneB_id=icmtid)
        sli_list.append(sli)

        raf1 = 'RAF1'
        raf1id = "NCBIGene:{}".format(self.entrez_dict.get(raf1))
        sli = self.get_sli(geneA_sym=nras, geneA_id=nrasid, geneB_sym=raf1, geneB_id=raf1id)
        sli_list.append(sli)
        sli = self.get_sli(geneA_sym=kras, geneA_id=krasid, geneB_sym=raf1, geneB_id=raf1id)
        sli_list.append(sli)

        shoc2 = 'SHOC2'
        shoc2id = "NCBIGene:{}".format(self.entrez_dict.get(shoc2))
        sli = self.get_sli(geneA_sym=nras, geneA_id=nrasid, geneB_sym=shoc2, geneB_id=shoc2id)
        sli_list.append(sli)
        sli = self.get_sli(geneA_sym=kras, geneA_id=krasid, geneB_sym=shoc2, geneB_id=shoc2id)
        sli_list.append(sli)

        prex1 = 'PREX1'
        prex1id = "NCBIGene:{}".format(self.entrez_dict.get(prex1))
        sli = self.get_sli(geneA_sym=nras, geneA_id=nrasid, geneB_sym=prex1, geneB_id=prex1id)
        sli_list.append(sli)
        sli = self.get_sli(geneA_sym=kras, geneA_id=krasid, geneB_sym=prex1, geneB_id=prex1id)
        sli_list.append(sli)

        return sli_list
