

class HgncEntry:
    def __init__(self, hgnc_id, symbol, entrez_id, ensembl_gene_id):
        self.hgnc_id = hgnc_id
        self.symbol = symbol
        self.entrez_id = entrez_id
        self.ensembl_gene_id = ensembl_gene_id

    def hgnc_id(self):
        return self.hgnc_id
    
    def symbol(self):
        return self.symbol
    
    def entrez_id(self):
        return self.entrez_id

    def ensembl_gene_id(self):
        return self.ensembl_gene_id