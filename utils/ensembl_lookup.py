#import ensembl_rest
from biomart import BiomartServer


import logging
import os
import gzip
from collections import defaultdict

import wget


class EnsemblLookup(object):
    """
    Ensembl Lookup service

    Uses Ensembl to build
    - Ensembl Gene ID to NCBI Entrez lookup
    - NCBI Entrez to Ensembl Gene ID lookup

    """
    LOOKUP_FILE = os.path.join(os.path.dirname(__file__), '..', 'lookup', 'gene2ensembl.gz')
    prot2gene_file = os.path.join(os.path.dirname(__file__), '..', 'lookup', 'prot2gene.gz')

    def __init__(self, filename=LOOKUP_FILE, backupfile = prot2gene_file,species_id=["9606"]):
        if not os.path.exists(filename):
            logging.warning("File {} does not exist.".format(filename))
            filename = self.download_file()

        self.filename = filename
        self.ncbi_lookup = defaultdict(str)
        self.symbol_lookup = defaultdict(str)
        self.reverse_ncbi_lookup = defaultdict(str)
        self.protein_lookup = defaultdict(str)
        self.backup_protein_lookup = defaultdict(str)

        with gzip.open(filename, 'rt') as f:
            for line in f:
                fields = line.split('\t')
                if fields[0].strip() in species_id:
                    ncbi_id = fields[1]
                    ensembl = fields[2]
                    prot_ID = fields[6].split(".")[0]
                    self.ncbi_lookup[ncbi_id] = ensembl
                    self.reverse_ncbi_lookup[ensembl] = ncbi_id
                    self.backup_protein_lookup[prot_ID] = ensembl

        with gzip.open(backupfile, 'rt') as f:
            next(f)  # skip header
            for line in f:
                fields = line.split('\t')
                if not fields[1] == "":
                    ensembl = fields[0]
                    prot_ID = fields[1]
                    self.protein_lookup[prot_ID] = ensembl
                    self.ncbi_lookup[ncbi_id] = ensembl
                    self.reverse_ncbi_lookup[ensembl] = ncbi_id

        with open(os.path.join(os.path.dirname(__file__), '..', 'lookup', 'lookup.txt'), 'rt') as f:
            next(f)  # skip header
            for line in f:
                fields = line.split('\t')
                symbol = fields[0]
                ncbi_id = fields[1]
                ensembl = fields[2]
                self.symbol_lookup[symbol] = ensembl
                self.ncbi_lookup[ncbi_id] = ensembl
                self.reverse_ncbi_lookup[ensembl] = ncbi_id


    def download_file(self,
                      urldir='ftp://ftp.ncbi.nih.gov/gene/DATA/',
                      filename='gene2ensembl.gz'):
        """
        Download the lookup file from NCBI FTP and save to its intended location.
        """
        url = os.path.join(urldir, filename)
        logging.info("Downloading file from {}".format(url))
        local_filename = wget.download(url, out=self.LOOKUP_FILE)
        return local_filename

    def map_ensembl_to_ncbiID(self, ensembl):
        """
        Map  a given Ensembl Gene ID to NCBI Entrez ID.
        """
        retval = None
        if ensembl in self.reverse_lookup:
            retval = self.reverse_lookup[ensembl]
        return retval

    def map_ncbiID_to_ensembl(self, ncbi_id):
        """
        Map a given Ensembl Gene ID to Ensembl Gene ID.
        """
        retval = None
        if ncbi_id in self.lookup:
            retval = self.lookup[ncbi_id]
        return retval

    def map_protID_to_ensembl(self, ncbi_id):
        """
        Map a given Ensembl Protein ID to Ensembl Gene ID.
        """
        retval = None
        if ncbi_id in self.lookup:
            retval = self.lookup[ncbi_id]
        return retval






#server = BiomartServer("http://www.ensembl.org/biomart")
#if server.is_alive:
#    print("ist am Leben")
#server.verbose = False
#ensembl = server.datasets['hsapiens_gene_ensembl']
#ensembl.show_filters()

#def ncbi2ensembl(list):
#    if not BiomartServer("http://www.ensembl.org/biomart").is_alive:
#        server = BiomartServer("http://www.ensembl.org/biomart")
#        ensembl = server.datasets['hsapiens_gene_ensembl']
#
#    res = []
#    response = ensembl.search({'filters': {'ensembl_peptide_id': list},
#                           'attributes': ['ensembl_gene_id']})
#    for line in response.iter_lines():
#        line = line.decode('utf-8')
#        res.append(line.split("\t")[0])
#    return res

#print(ncbi2ensembl(['79577', '5981', '641', '990', '7465', '5884', '51343', '1111']))
#['ENSG00000149554', 'ENSG00000105325', 'ENSG00000276618', 'ENSG00000152942', 'ENSG00000035928', 'ENSG00000197299', 'ENSG00000166483', 'ENSG00000134371', 'ENSG00000094804']
#  CHEK1              FZR1               RAD17











