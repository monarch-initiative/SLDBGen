from __future__ import print_function
import gzip
from mailbox import FormatError
import pandas as pd
import os.path
from collections import defaultdict
from utils.lookup import Lookup
import json


import logging.handlers

handler = logging.handlers.WatchedFileHandler(
    os.environ.get("LOGFILE", "parser.log"))
formatter = logging.Formatter('%(asctime)s - %(levelname)s -%(filename)s:%(lineno)d - %(message)s')
handler.setFormatter(formatter)
log = logging.getLogger()
log.setLevel(os.environ.get("LOGLEVEL", "INFO"))
log.addHandler(handler)




class Parser:
    """
    A class to parse the graph.
    """

    def __init__(self, data_dir="data", data_lookup_dir ="lookup", params={}):

        if os.path.exists(data_dir) and not os.path.isdir(data_dir):
            raise Exception("`data_dir` must point to a directory")
        elif not os.path.exists(data_dir):
            raise Exception("`data_dir` does not exist. Make it and download the required files to run this script.")
        self.data_dir = data_dir
        self.data_lookup_dir = data_lookup_dir
        log.info("Initializing Parser with data directory: {}".format(self.data_dir))
        # Set up values for parameters
        default_string_ppi_threshold = 700
        # the following use defaults if params does not contain the item
        self.string_interaction_threshold_ = params.get('string_threshold', default_string_ppi_threshold)

        # Same for files. Note that the user can override the defaults by providing arguments for the paths
        # to these files.
        # latest files
        default_gene2ensembl_path = os.path.join(data_dir, "gene2ensembl.gz")
        default_string_path = os.path.join(data_dir, "STRING_9606.protein.links.v11.0.txt.gz")
        default_prot2gene_path =  os.path.join(data_dir, "STRING_9606.protein.links.v11.0.txt.gz")
        # Dictionary of main results to print out for user
        self.summary = defaultdict(str)
        self.gene2ensembl_path = params.get('gene2ensembl_path', default_gene2ensembl_path)
        self.string_ppi_path = params.get('string_path', default_string_path)
        self.prot2gene_path = params.get('prot2gene_path', default_prot2gene_path)


        # check whether the files actually exist. If not, raise an Error
        self.qc_input()

        log.info("[INFO] Parameters:")
        log.info("\tString interaction threshold: {}".format(self.string_interaction_threshold_))

    def qc_input(self):
        if not os.path.isfile(self.gene2ensembl_path):
            raise Exception("Could not find gene2ensembl file at {}".format(self.gene2ensembl_path))
        else:
            log.info("gene2ensembl: {}".format(self.gene2ensembl_path))
        if not os.path.isfile(self.string_ppi_path):
            raise Exception("Could not find STRING PPI file at {}".format(self.string_ppi_path))
        else:
            log.info("STRING PPI file: {}".format(self.string_ppi_path))





    def load_lookup(self):
        protein2ensembl = defaultdict(str)
        with gzip.open(self.prot2gene_path, 'rt') as f:
            next(f)  # skip header
            for line in f:
                fields = line.split('\t')
                if fields[0].strip() in ["9606"]:
                    ensembl = fields[2]
                    if fields[6].startswith("E") and fields[6] is not None:
                        prot_ID = fields[6].split(".")[0]
                        protein2ensembl[prot_ID] = ensembl

        return protein2ensembl

    def load_backup(self):
        protein2ensembl = defaultdict(str)
        with gzip.open(self.prot2gene_path, 'rt') as f:
            next(f)  # skip header
            for line in f:
              fields = line.split('\t')
              ensembl = fields[0]
              if fields[1].startswith("E") and fields[1] is not None:
                prot_ID = fields[1]
                protein2ensembl[prot_ID] = ensembl
        return protein2ensembl
