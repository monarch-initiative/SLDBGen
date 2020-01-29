import os
import json
import gzip
from collections import defaultdict

import wget


class EntrezLookup(object):
    LOOKUP_FILE = os.path.join(os.path.dirname(__file__), '..', 'lookup', 'Homo_sapiens.gene_info.gz')

    def __init__(self, filename=LOOKUP_FILE):
        if not os.path.exists(filename):
            print("File {} does not exist.".format(filename))
            filename = self.download_file()

        self.filename = filename
        self.lookup = defaultdict(str)
        self.reverse_lookup = defaultdict(str)
        self.dbxrefs = defaultdict(str)

        with gzip.open(filename, 'rt') as f:
            for line in f:
                if line.startswith("9606"):
                    fields = line.split('\t')
                    gene_id = fields[1]
                    symbol = fields[2]
                    dbxrefs = fields[5]
                    self.lookup[gene_id] = symbol
                    self.reverse_lookup[symbol] = gene_id
                    self.dbxrefs[gene_id] = dbxrefs.split('|')

    def download_file(self):
        urldir = 'ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/'
        filename = 'Homo_sapiens.gene_info.gz'
        if not os.path.exists(self.LOOKUP_FILE):
            url = os.path.join(urldir, filename)
            print("Downloading file from {}".format(url))
            local_filename = wget.download(url, out=self.LOOKUP_FILE)
        return local_filename

    def map_symbol_to_identifier(self, symbol):
        retval = None
        if symbol in self.reverse_lookup:
            retval = self.reverse_lookup[symbol]
        return retval

    def map_identifier_to_symbol(self, identifier):
        retval = None
        if identifier in self.lookup:
            retval = self.lookup[identifier]
        return retval

    def map_symbol_to_specific_identifier(self, symbol, prefix):
        identifier = self.map_symbol_to_identifier(symbol)
        dbxrefs = self.dbxrefs[identifier]
        filtered = list(filter(lambda x: (prefix in x), dbxrefs))
        return filtered

