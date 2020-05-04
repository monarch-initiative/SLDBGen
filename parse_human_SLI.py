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
yeastSymbol2entrezID = EntrezLookup(filename="lookup/Saccharomyces_cerevisiae.gene_info.gz",
                                    species_id=["4932", "559292"]
                                    ).reverse_lookup


# Do the yeast somewhere else
# boone_sli_list = idg2sl.parse_costanzo_boone_2016_NxN_data(symbol2id=yeastSymbol2entrezID)

# luo2008 = idg2sl.parse_luo2009_supplemental_file_S3('data/luo2009.tsv', humanSymbol2entrezID)
# bommi2008 = idg2sl.parse_bommi_reddi_2008('data/bommi-reddy-2008.tsv', humanSymbol2entrezID)
# turner_list = idg2sl.parse_turner_2008('data/turner-PARP1-2008.tsv', humanSymbol2entrezID)

# steckel2012 = idg2sl.parse_steckel_2012('data/steckel-2012-KRAS.tsv', humanSymbol2entrezID)
# lord2008 = idg2sl.parse_lord_2008('data/lord-PARP1-2008.tsv', humanSymbol2entrezID)
# toyoshima2008 = idg2sl.parse_toyoshima_2008('data/toyoshima-MYC-2008.tsv', humanSymbol2entrezID)
# shen2015 = idg2sl.parse_Shen2015('data/Shen_2015.txt', humanSymbol2entrezID)

# srivas2016 = idg2sl.parse_srivas_2016('data/Srivas_2016.txt', humanSymbol2entrezID)

han2017 = idg2sl.parse_han_2017('data/Han2017_supplemental_table_4.txt', humanSymbol2entrezID)

# sli_lists = [luo2008, bommi2008, turner_list, steckel2012, lord2008, toyoshima2008, shen2015, srivas2016]
sli_lists = [han2017]

n = 0
n_SL = 0
for sli_list in sli_lists:
    for sli in sli_list:
        print(sli.get_tsv_line())
        n += 1
        if sli.get_SL():
            n_SL += 1

print("We got %d interactions including %d synthetic lethal interactions" % (n, n_SL))


