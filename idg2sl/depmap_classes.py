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