class Gene:
    """
    Instances of this class represent a gene from the depmap data.
    Each gene has a a set of cells it is effecitve in.
    """

    def __init__(self,
                 gene_id="",
                 effective_cells=None,
                 proportion_effective_cells=None):
        self.gene_id = gene_id
        self.effective_cells = effective_cells
        self.proportion_effective_cells = proportion_effective_cells

    def get_id(self):
        return self.gene_id

    def get_effective_cells(self):
        return self.effective_cells

    def get_proportion(self):
        return self.proportion_effective_cells

    def get_number_effective_cells(self):
        return len(self.effective_cells)


class SyntheticLethalInteraction:
    """
    Instances of this class represent a single synthetic lethality (SL) interaction
    """

    def __init__(self,
                 gene_A_id="",
                 gene_B_id="",
                 SL=None):
        self.gene_A_id = gene_A_id
        self.gene_B_id = gene_B_id
        self.SL = SL

    def get_gene_A_id(self):
        return self.gene_A_id

    def get_gene_B_id(self):
        return self.gene_B_id

    def get_SL(self):
        return self.SL
