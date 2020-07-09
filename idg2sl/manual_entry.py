

breastCarcinoma = "Breast Carcinoma"
ncitBreastCarcinoma = "NCIT:C4872"
growthInhibitionAssay = 'growth.inhibition.assay'


## In a number of papers, only one or a handful of synthetic lethal interactions are described.
## These are very valuable. It is easiest to enter this information by hand.

class ManualEntry:
    def __init__(self):
        self.entries = []

    def _get_Reid2016(self):
        """
         PLK1 NCBI Gene id 5347 and CKS1B NCBI Gene id 1163 
         """
        sli = SyntheticLethalInteraction(gene_A_symbol='PLK1',
                                     gene_A_id=5347,
                                     gene_B_symbol = "CKS1B",
                                     gene_B_id = 1163,
                                     gene_A_pert='shRNA',
                                     gene_B_pert='increased.expression',
                                     effect_type='growth.inhibition',
                                     effect_size = 'n/a',
                                     cell_line='multiple.breast.cancer.cell.lines',
                                     cellosaurus_id='n/a',
                                     cancer_type=breastCarcinoma,
                                     ncit_id=ncitBreastCarcinoma,
                                     assay=growthInhibitionAssay,
                                     pmid='PMID:27558135',
                                     SL=True)

    def get_entries(self):
        self.entries.append(self._get_Reid2016())
    

