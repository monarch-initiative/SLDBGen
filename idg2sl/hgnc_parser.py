import csv
from collections import defaultdict
from .hgnc_entry import HgncEntry

class HgncParser:
    def __init__(self, fname):
        self.fname = fname
        self.d = defaultdict(HgncEntry)
        with open(fname, 'r') as csvfile: 
            csvreader = csv.DictReader(csvfile,delimiter='\t') 
            for row in csvreader: 
                hgnc_id = row['hgnc_id'] 
                symbol  = row['symbol']
                entrez_id = row['entrez_id']
                ensembl_gene_id = row['ensembl_gene_id']
                self.d[symbol] = HgncEntry(hgnc_id, symbol, entrez_id, ensembl_gene_id)

    def get_dictionary(self):
        return self.d
