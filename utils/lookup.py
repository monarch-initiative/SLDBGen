import logging
import os
import json
import gzip
from collections import defaultdict

import wget


class Lookup(object):
    """
    Lookup service

    Builds
    - NCBI Gene ID to symbol lookup
    - NCBI Gene ID to Ensembl Gene ID lookup
    - Ensembl Gene ID to dbxref lookup
    - Ensembl Gene ID to Symbol
    - Symbol to EnsemblGene ID lookup
    - Symbol to NCBI Gene ID lookup

    """
    LOOKUP_FILE = os.path.join(os.path.dirname(__file__), '..', 'lookup', 'Homo_sapiens.gene_info.gz')

    def __init__(self, filename=LOOKUP_FILE, species_id=["9606"]):
        if not os.path.exists(filename):
            logging.warning("File {} does not exist.".format(filename))
            filename = self.download_file()

        self.filename = filename

        self.ncbi2symbol = defaultdict(str)
        self.ncbi2ensembl = defaultdict(str)
        self.symbol2ncbi = defaultdict(str)
        self.symbol2ensembl = defaultdict(str)
        self.ensembl2symbol = defaultdict(str)
        self.ensembl2ncbi = defaultdict(str)
        self.protein2ensembl = defaultdict(str)
        self.dbxrefs = defaultdict(str)
        #dicts = [self.symbol2ensembl, self.symbol2ncbi,
        #     self.ensembl2symbol, self.ensembl2ncbi,
        #     self.ncbi2ensembl, self.ncbi2symbol,
        #     self.protein2ensembl, self.dbxrefs]

        ### Homo_sapiens.gene_info.gz
        with gzip.open(filename, 'rt') as f:
            for line in f:
                fields = line.split('\t')
                if fields[0].strip() in species_id:
                    ncbi = fields[1]
                    symbol = fields[2]
                    ncbi_synonyms = fields[4].split("|")
                    dbxrefs = fields[5]
                    ensembl = dbxrefs.split(":")[-1]
                    self.ncbi2symbol[ncbi] = symbol
                    self.dbxrefs[ncbi] = dbxrefs.split('|')
                    self.symbol2ncbi[symbol] = ncbi
                    for syn in ncbi_synonyms:
                        self.symbol2ncbi[syn] = ncbi
                    if ensembl.startswith("ENS") and ensembl is not None:
                        self.ncbi2ensembl[ncbi] = ensembl
                        self.ensembl2ncbi[ensembl] = ncbi
                        self.ensembl2symbol[ensembl] = symbol
                        self.symbol2ensembl[symbol] = ensembl
                        for syn in ncbi_synonyms:
                            self.symbol2ensembl[syn] = ensembl

        ### lookup.txt
        with open(os.path.join(os.path.dirname(__file__), '..', 'lookup', 'lookup.txt'), 'rt') as f:
            next(f)  # skip header
            for line in f:
                fields = line.split('\t')
                symbol = fields[0]
                if fields[1] is not None:
                    ncbi = fields[1]
                else:
                    ncbi = "n/a"
                if fields[2] is not None:
                    ensembl = fields[2]
                else:
                    ensembl = "n/a"

                if ncbi not in self.ncbi2ensembl.keys() and ensembl is not "n/a":
                    self.ncbi2ensembl[ncbi] = ensembl
                if ncbi not in self.ncbi2symbol.keys() and symbol is not "n/a":
                    self.ncbi2symbol[ncbi] = symbol
                if ensembl not in self.ensembl2ncbi.keys() and ncbi is not "n/a":
                    self.ensembl2ncbi[ensembl] = ncbi
                if ensembl not in self.ensembl2symbol.keys() and symbol is not "n/a":
                    self.ensembl2symbol[ensembl] = symbol
                if symbol not in self.symbol2ensembl.keys() and ensembl is not "n/a":
                    self.symbol2ensembl[symbol] = ensembl
                if symbol not in self.symbol2ncbi.keys() and ncbi is not "n/a":
                    self.symbol2ncbi[symbol] = ncbi

        with gzip.open(os.path.join(os.path.dirname(__file__), '..', 'lookup', 'gene2ensembl.gz'), 'rt') as f:
            next(f)  # skip header
            for line in f:
                fields = line.split('\t')
                if fields[0].strip() in species_id:
                    ncbi = fields[1]
                    ensembl = fields[2]
                    if fields[6].startswith("E") and fields[6] is not None:
                        prot_ID = fields[6].split(".")[0]
                        #print(prot_ID)
                        if prot_ID not in self.protein2ensembl.keys():
                            self.protein2ensembl[prot_ID] = ensembl

                    if ncbi not in self.ncbi2ensembl.keys():
                        self.ncbi2ensembl[ncbi] = ensembl
                    if ensembl not in self.ensembl2ncbi.keys():
                        self.ensembl2ncbi[ensembl] = ncbi

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
        if symbol in self.symbol2ncbi:
            retval = self.symbol2ncbi[symbol]
        return retval

    def map_identifier_to_symbol(self, identifier):
        """
        Map a given NCBI Gene ID to symbol.
        """
        retval = None
        if identifier in self.ncbi2symbol:
            retval = self.ncbi2symbol[identifier]
        return retval

    def map_symbol_to_specific_identifier(self, symbol, prefix):
        """
        Map a given symbol to an dbxref identifier that has the given prefix.
        """
        identifier = self.map_symbol_to_identifier(symbol)
        dbxrefs = self.dbxrefs[identifier]
        filtered = list(filter(lambda x: (prefix in x), dbxrefs))
        return filtered

