import logging
import os
import json
import gzip
from collections import defaultdict

import wget


class EntrezLookup(object):
    """
    Entrez Lookup service

    Uses NCBI Entrez to build
    - NCBI Gene ID to symbol lookup
    - symbol to NCBI Gene ID lookup
    - NCBI Gene ID to dbxref lookup

    """
    LOOKUP_FILE = os.path.join(os.path.dirname(__file__), '..', 'lookup', 'Homo_sapiens.gene_info.gz')

    def __init__(self, filename=LOOKUP_FILE, species_id=["9606"]):
        if not os.path.exists(filename):
            logging.warning("File {} does not exist.".format(filename))
            filename = self.download_file()

        self.filename = filename
        self.lookup = defaultdict(str)
        self.reverse_lookup = defaultdict(str)
        self.dbxrefs = defaultdict(str)

        with gzip.open(filename, 'rt') as f:
            for line in f:
                fields = line.split('\t')
                if fields[0].strip() in species_id:
                    gene_id = fields[1]
                    symbol = fields[2]
                    dbxrefs = fields[5]
                    self.lookup[gene_id] = symbol
                    self.reverse_lookup[symbol] = gene_id
                    self.dbxrefs[gene_id] = dbxrefs.split('|')

    def download_file(self,
                      urldir='ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/',
                      filename='Homo_sapiens.gene_info.gz'):
        """
        Download the lookup file from NCBI FTP and save to its intended location.
        """
        url = os.path.join(urldir, filename)
        logging.info("Downloading file from {}".format(url))
        local_filename = wget.download(url, out=self.LOOKUP_FILE)
        return local_filename

    def map_symbol_to_identifier(self, symbol):
        """
        Map  a given symbol to NCBI Gene ID.
        """
        retval = None
        if symbol in self.reverse_lookup:
            retval = self.reverse_lookup[symbol]
        return retval

    def map_identifier_to_symbol(self, identifier):
        """
        Map a given NCBI Gene ID to symbol.
        """
        retval = None
        if identifier in self.lookup:
            retval = self.lookup[identifier]
        return retval

    def map_symbol_to_specific_identifier(self, symbol, prefix):
        """
        Map a given symbol to an dbxref identifier that has the given prefix.
        """
        identifier = self.map_symbol_to_identifier(symbol)
        dbxrefs = self.dbxrefs[identifier]
        filtered = list(filter(lambda x: (prefix in x), dbxrefs))
        return filtered

