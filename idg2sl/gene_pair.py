
class GenePair:
    """
    This class acts as a key in dictionaries with the key being made of gene A and gene B
    because in some of our datasets, both gene A and gene B can vary
    """
    def __init__(self, geneA, geneB):
        self.gene_A = geneA
        self.gene_B = geneB

    def __hash__(self):
        return hash((self.gene_A, self.gene_B))

    def __eq__(self, other):
        return (self.gene_A, self.gene_B) == (other.gene_A, other.gene_B)

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        return not (self == other)
