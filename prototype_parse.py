import gzip
import os.path
import wget

from collections import defaultdict
from utils.entrez_lookup import EntrezLookup
import idg2sl


def get_entrez_gene_map():
    """
    Download the file  ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
    """
    urldir = 'ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/'
    local_filename = 'Homo_sapiens.gene_info.gz'
    symbol2entrezID = defaultdict(str)
    if not os.path.exists(local_filename):
        url = os.path.join(urldir, local_filename)
        local_filename = wget.download(url)
    with gzip.open(local_filename, 'rt') as f:
        for line in f:
            if not line.startswith("9606"):
                continue  # non human homo sapiens
            fields = line.split('\t')
            entrez = fields[1]
            symbol = fields[2]
            symbol2entrezID[symbol] = entrez
    return symbol2entrezID


humanSymbol2entrezID = EntrezLookup().reverse_lookup
yeastSymbol2entrezID = EntrezLookup(filename=
                                 "lookup/Saccharomyces_cerevisiae.gene_info.gz",
                                    species_id=["4932", "559292"]
                                    ).reverse_lookup

luo_list = idg2sl.parse_luo2009_supplemental_file_S3('data/luo2009.tsv', humanSymbol2entrezID)
bommi_list = idg2sl.parse_bommi_reddi_2008('data/bommi-reddy-2008.tsv', humanSymbol2entrezID)
turner_list = idg2sl.parse_turner_2008('data/turner-PARP1-2008.tsv', humanSymbol2entrezID)
boone_sli_list = idg2sl.parse_costanzo_boone_2016_NxN_data(symbol2id=yeastSymbol2entrezID)

sli_lists = [luo_list, bommi_list, turner_list, boone_sli_list]

for sli_list in sli_lists:
    for sli in sli_list:
        print(sli.get_tsv_line())


