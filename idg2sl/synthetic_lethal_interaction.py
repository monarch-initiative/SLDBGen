
class SyntheticLethalInteraction:
    """
    Instances of this class represent a single synthetic lethality (SL) interaction
    together with data about the experiment that was used to show the SL. A series
    of ingest scripts will transform the data of each of the ca. 30 published
    experiments to this common datastructure
    """

    def __init__(self, gene_A_symbol=None,
                 gene_A_id=None,
                 gene_B_symbol=None,
                 gene_B_id=None,
                 gene_A_pert=None,
                 gene_B_pert=None,
                 effect_type=None,
                 effect_size=None,
                 assay=None,
                 pmid=None,
                 SL=None):
        if gene_A_symbol is None:
            raise ValueError("Need to pass gene A")
        if gene_A_id is None:
            raise ValueError("Need to pass gene A id")
        if gene_B_symbol is None:
            raise ValueError("Need to pass gene B")
        if gene_B_id is None:
            raise ValueError("Need to pass gene B id")
        if gene_A_pert is None:
            raise ValueError("Need to pass gene A perturbation")
        if gene_B_pert is None:
            raise ValueError("Need to pass gene B perturbation")
        if assay is None:
            raise ValueError("Need to pass assay")
        if pmid is None:
            raise ValueError("Need to pass pmid")
        if SL is None:
            raise ValueError("Need to pass True or False for SL")
        self.gene_A_symbol = gene_A_symbol
        self.gene_A_id = gene_A_id
        self.gene_B_symbol = gene_B_symbol
        self.gene_B_id = gene_B_id
        self.gene_A_pert = gene_A_pert
        self.gene_B_pert = gene_B_pert
        self.assay = assay
        self.pmid = pmid
        if effect_type is None or effect_size is None:
            self.effect_type = ""
            self.effect_size = ""
        else:
            self.effect_type = effect_type
            self.effect_size = effect_size
        self.SL = SL # True: synthetic lethal, False: negative control

    def get_tsv_line(self):
        if self.SL:
            sl = 'T'
        else:
            sl = 'F'
        lst = [self.gene_A_symbol,
               str(self.gene_A_id),
               self.gene_B_symbol,
               str(self.gene_B_id),
               self.gene_A_pert,
               self.gene_B_pert,
               self.effect_type,
               str(self.effect_size),
               self.assay,
               self.pmid,
               sl]
        return "\t".join(lst)
