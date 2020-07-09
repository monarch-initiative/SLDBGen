from collections import defaultdict
from idg2sl.gene_pair import GenePair
from idg2sl.synthetic_lethal_interaction import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser

class Turner2008Parser(SL_DatasetParser):
    """
    Parse data from  Turner NC, et al., A synthetic lethal siRNA screen identifying genes mediating
    sensitivity to a PARP inhibitor. EMBO J. 2008 May 7;27(9):1368-77. PubMed PMID: 18388863
    PARP1 was inhibited by the PARP inhibitor KU0058948, and a short interfering RNA library targeting
    779 human protein kinase and kinase assocaited genes was applied.
    :param path:
    """
    def __init__(self, symbol2entrezID, fname, hasHeader=True):
        fname = fname
        pmid = 'PMID:18388863'

        super().__init__(fname=fname, pmid=pmid, symbol2entrezID = symbol2entrezID, hasHeader=hasHeader)

    def parse(self):
        print("Number of rows ", len(self.rows))
        parp1_symbol = 'PARP1'
        parp1_id = 'NCBIGene:142'
        parp1_perturbation = 'drug'
        gene2_perturbation = 'siRNA'
        assays = ['competitive hybridization', 'multicolor competition assay']
        assay_string = ";".join(assays)
        effect_type = 'stddev'
        cell_line = 'CAL-51'
        cellosaurus = 'CVCL_1110'
        cancer = "Breast Carcinoma"
        ncit = "NCIT:C4872"
        sli_dict = defaultdict(list)
        for row in self.rows:
            if len(row) < 3:
                raise ValueError("Bad row for Turner et al, only %d fields found" % len(row))
            geneB_sym = row[0]
            if geneB_sym in self.symbol2entrezID:
                geneB_id = "NCBIGene:{}".format(self.symbol2entrezID.get(geneB_sym))
            else:
                geneB_id = "n/a"
            zscore = float(row[1])
            percent_inhib = row[2]
            if zscore <= -3.0:
                SL = True
            else:
                SL = False
            sli = SyntheticLethalInteraction(gene_A_symbol=parp1_symbol,
                                             gene_A_id=parp1_id,
                                             gene_B_symbol=geneB_sym,
                                             gene_B_id=geneB_id,
                                             gene_A_pert=parp1_perturbation,
                                             gene_B_pert=gene2_perturbation,
                                             effect_type=effect_type,
                                             effect_size=zscore,
                                             cell_line=cell_line,
                                             cellosaurus_id=cellosaurus,
                                             cancer_type=cancer,
                                             ncit_id=ncit,
                                             assay=assay_string,
                                             pmid=self.pmid,
                                             SL=SL)
            gene_pair = GenePair(parp1_symbol, geneB_sym)
            sli_dict[gene_pair].append(sli)
            sli_list = self._mark_maximum_entries(sli_dict)
        return sli_list
