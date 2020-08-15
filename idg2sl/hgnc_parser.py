import csv
from collections import defaultdict

class HgncParser:
    def __init__(self, fname):
        self.fname = fname
        self.symbol2entrez = defaultdict(str)
        self.symbol2ensembl = defaultdict(str)
        self.synonym2symbol = defaultdict(str)
        with open(fname, 'r') as csvfile: 
            csvreader = csv.DictReader(csvfile, delimiter='\t')
            for row in csvreader: 
                # hgnc_id = row['hgnc_id']
                symbol = row['symbol']
                entrez_id = row['entrez_id']
                ensembl_gene_id = row['ensembl_gene_id']
                prevsymstr = row['prev_symbol'].replace('"', "")
                prev_symbols = prevsymstr.split('|')
                for ps in prev_symbols:
                    if ps in self.synonym2symbol:
                        self.synonym2symbol[ps] = "MULTIPLE"
                    else:
                        self.synonym2symbol[ps] = symbol
                alias_symbol_str = row['alias_symbol'].replace('"', "")
                alias_symbols = alias_symbol_str.split('|')
                for ps in alias_symbols:
                    if ps in self.synonym2symbol:
                        self.synonym2symbol[ps] = "MULTIPLE"
                    else:
                        self.synonym2symbol[ps] = symbol
                self.symbol2entrez[symbol] = entrez_id
                self.symbol2ensembl[symbol] = ensembl_gene_id

    def get_entrez_dictionary(self):
        return self.symbol2entrez

    def get_ensembl_dictionary(self):
        return self.symbol2ensembl

    def get_synonym_dictionary(self):
        return self.synonym2symbol