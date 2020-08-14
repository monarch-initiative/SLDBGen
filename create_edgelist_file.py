# Can be run from ipython as %run  run_hnn2v_analysis.py
import sys
import os
from parser_all_files import Parser


import logging.handlers
handler = logging.handlers.WatchedFileHandler(
    os.environ.get("LOGFILE", "parser.log"))
formatter = logging.Formatter('%(asctime)s - %(levelname)s -%(filename)s:%(lineno)d - %(message)s')
handler.setFormatter(formatter)
log = logging.getLogger()
log.setLevel(os.environ.get("LOGLEVEL", "INFO"))
log.addHandler(handler)



data_download_dir = "data"
# 1) Check the required files are present in data directory
def check_existence(url, local_filepath):
    if not os.path.isfile(local_filepath):
        log.error("could file not find at {}".format(local_filepath))
        log.error("[FATAL] Need to download file from {} and place the file in the directory.".format(url))
        sys.exit(1)
    else:
        log.info("Using file {}".format(local_filepath))

string_protein_links_url = "https://stringdb-static.org/download/protein.links.v11.0/9606.protein.links.v11.0.txt.gz"
string_local_filename = "9606.protein.actions.v11.0.txt.gz"
string_path = os.path.join(data_download_dir, string_local_filename)
check_existence(string_protein_links_url, string_path)

depmap_url = 'https://ndownloader.figshare.com/files/14221385'
filename = 'data/gene_effect_corrected.csv'  # path to the file

data_lookup_dir = "lookup"

gene2ensembl_url = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz"
gene2ensembl_local_file = "gene2ensembl.gz"
gene2ensembl_path = os.path.join(data_download_dir, gene2ensembl_local_file)
check_existence(gene2ensembl_url, gene2ensembl_path)

prot2gene_url = "" #TODO: add the url for this file
prot2gene_file = "prot2gene.gz"
prot2gene_path = os.path.join(data_lookup_dir, prot2gene_file)
check_existence(prot2gene_url, prot2gene_path)



parameters = {'string_threshold': 700,
              'prot2gene_path': prot2gene_path,
              'string_path': string_path,
              'gene2ensemble_path':gene2ensembl_path
              }

parser = Parser(params=parameters)
parser.parse()
parser.print_summary()
outfilename_training = "edgelist.txt"
parser.gene_gene_intercations(outfilename_training)
