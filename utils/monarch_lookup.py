import os
import json
import gzip


class MonarchLookup(object):
    LOOKUP_FILE = os.path.join(os.path.dirname(__file__), '..', 'lookup', 'monarch_human_gene_lookup.json.gz')

    def __init__(self, filename=LOOKUP_FILE):
        if not os.path.exists(filename):
            raise TypeError("File {} does not exist.".format(filename))

        self.filename = filename
        self.lookup = {}
        self.reverse_lookup = {}
        self.dbxrefs = {}

        with gzip.open(filename, "rb") as f:
            file_json = json.loads(f.read().decode("utf-8"))

        docs = file_json['response']['docs']
        for d in docs:
            self.lookup[d['id']] = d['label'][0]
            self.reverse_lookup[d['label'][0]] = d['id']
            if 'equivalent_curie' in d:
                self.dbxrefs[d['id']] = d['equivalent_curie']

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
        equivalent_curies = self.dbxrefs[identifier]
        filtered = list(filter(lambda x: (prefix in x), equivalent_curies))
        return filtered

