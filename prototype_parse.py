import csv
import gzip
import logging
import os.path
from collections import defaultdict

import requests
import zipfile
import tempfile

from utils.entrez_lookup import EntrezLookup
from clint.textui import progress


class SyntheticLethalInteraction:
    """
    Instances of this class represent a single synthetic lethality (SL) interaction
    together with data about the experiment that was used to show the SL. A series
    of ingest scripts will transform the data of each of the ca. 30 published
    experiments to this common datastructure
    """

    def __init__(self, gene_A_symbol=None,
                 gene_A_id=None,
                 gene_B_symbol=None,
                 gene_B_id=None,
                 gene_A_pert=None,
                 gene_B_pert=None,
                 effect_type=None,
                 effect_size=None,
                 assay=None,
                 pmid=None):
        if gene_A_symbol is None:
            raise ValueError("Need to pass gene A")
        if gene_A_id is None:
            raise ValueError("Need to pass gene A id")
        if gene_B_symbol is None:
            raise ValueError("Need to pass gene B")
        if gene_B_id is None:
            raise ValueError("Need to pass gene B id")
        if gene_A_pert is None:
            raise ValueError("Need to pass gene A perturbation")
        if gene_B_pert is None:
            raise ValueError("Need to pass gene B perturbation")
        if assay is None:
            raise ValueError("Need to pass assay")
        if pmid is None:
            raise ValueError("Need to pass pmid")
        self.gene_A_symbol = gene_A_symbol
        self.gene_A_id = gene_A_id
        self.gene_B_symbol = gene_B_symbol
        self.gene_B_id = gene_B_id
        self.gene_A_pert = gene_A_pert
        self.gene_B_pert = gene_B_pert
        self.assay = assay
        self.pmid = pmid
        if effect_type is None or effect_size is None:
            self.effect_type = ""
            self.effect_size = ""
        else:
            self.effect_type = effect_type
            self.effect_size = effect_size

    def get_tsv_line(self):
        lst = [self.gene_A_symbol,
               str(self.gene_A_id),
               self.gene_B_symbol,
               str(self.gene_B_id),
               self.gene_A_pert,
               self.gene_B_pert,
               self.effect_type,
               str(self.effect_size),
               self.assay,
               self.pmid]
        return "\t".join(lst)


def parse_luo2009_supplemental_file_S3(path, symbol2entrezID):
    """
    Luo J, et al A genome-wide RNAi screen identifies multiple synthetic lethal
    interactions with the Ras oncogene. Cell. 2009 May 29;137(5):835-48.
    PubMed PMID: 19490893
    Supplemental Table S3 was copied from the original PDF file and stored to
    a file 'luo2009.tsv' in the data directory. The experiment used cell lines
    with activating (oncogenic) mutations in KRAS and used a first screen with
    relative abundance of shRNAs to identify 368 genes using stringent criteria.
    They tested 320 candidates from the frist screen and found 83 shRNAs targeting
    77 genes to preferentially decreased the viability of the KRAS mutant cells
    compared to WT cells, thus indicating SL. They screened 68 of these candidates
    in a second line and could confirm 50 of them (73.5%), indicating that the
    majority of the candidates were likely to be true positives.
    """
    # because of the experiment, geneA is always KRAS. GeneB is in Table S3
    kras_symbol = 'KRAS'
    kras_id = 'NCBIGene:3845'
    kras_perturbation = 'oncogenic_mutation'
    gene2_perturbation = 'shRNA'
    pmid = 'PMID:19490893'
    assays = ['competitive hybridization', 'multicolor competition assay']
    assay_string = ";".join(assays)
    effect_type = 'stddev'

    if not os.path.exists(path):
        raise ValueError("Must path a valid path for Luo et al 2009")
    sli_dict = defaultdict(SyntheticLethalInteraction)

    with open(path) as f:
        next(f)  # skip header
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 8:
                logging.error("Only got %d fields but was expecting 8" % len(fields))
                i = 0
                for f in fields:
                    print("%d) %s" % (i, f))

                raise ValueError("Malformed line, must have 8 tab-separated fields")

            geneB_sym = fields[0]
            geneB_refSeq = fields[1]
            if geneB_sym in symbol2entrezID:
                geneB_id = "NCBIGene:{}".format(symbol2entrezID.get(geneB_sym))
            else:
                geneB_id = 'n/a'
            stddev = float(fields[5])
            sli = SyntheticLethalInteraction(gene_A_symbol=kras_symbol,
                                             gene_A_id=kras_id,
                                             gene_B_symbol=geneB_sym,
                                             gene_B_id=geneB_id,
                                             gene_A_pert=kras_perturbation,
                                             gene_B_pert=gene2_perturbation,
                                             effect_type=effect_type,
                                             effect_size=stddev,
                                             assay=assay_string,
                                             pmid=pmid)
            if geneB_sym in sli_dict:
                # get the entry with the strongest effect size
                sli_b = sli_dict.get(geneB_sym)
                if abs(stddev) > abs(sli.effect_size):
                    sli_dict[geneB_sym] = sli
            else:
                # first entry for geneB
                sli_dict[geneB_sym] = sli
    return sli_dict.values()


def parse_costanzo_boone_2016_NxN_data(symbol2id,
                                       force_download=False,
                                       pvalue_cutoff=0.05
                                       ) -> defaultdict:
    """
    Costanzo et al. A global genetic interaction network maps a wiring diagram of
    cellular function. Science. 23 Sep 2016: Vol. 353. Issue 6306.
    DOI: 10.1126/science.aaf1420

    This paper describes the results of massive S. cerevisiae SGA experiment.
    The method parses data describing the effect of knocking out all pairwise
    combinations of non-essential genes (SGA_NxN.txt).

    `This method parses a data file (SGA_NxN.txt) extracted from this zip file:
    <http://boonelab.ccbr.utoronto.ca/supplement/costanzo2016/data_files/
    Data%20File%20S1_Raw%20genetic%20interaction%20datasets:%20Pair-wise%20
    interaction%20format.zip>`_

    This method saves SLIs with pvalue < pvalue_cutoff. Downstream user can filter
    further using epsilon score. See methods for details:
        "The interaction data ... should be filtered prior to use. We suggest three
        different thresholds [lenient (P < 0.05), intermediate (P < 0.05 and e < 0.08),
        and stringent confidence (P <0.05 and e > 0.16 or e < -0.12)]"

    `Detailed description of Supplemental Data is here:
    <http://boonelab.ccbr.utoronto.ca/supplement/costanzo2016/>`_

    :param symbol2id: dict produced by `EntrezLookup` to look up IDs from symbols
    :param force_download: force dl of zip file and rewriting interaction data [false]
    :param pvalue_cutoff: only SLIs with pvalues less than this will be output [0.05]
    :return: defaultdict with SL interactions
    """
    interaction_data_file = "Data File S1. Raw genetic interaction datasets: " \
                            "Pair-wise interaction format/" \
                            "SGA_NxN.txt"
    local_data_file = 'data/SGA_NxN.txt.gz'  # zip it up after download, it's big
    zip_file_url = 'http://boonelab.ccbr.utoronto.ca/supplement/costanzo2016/' \
                   'data_files/' \
                   'Data%20File%20S1_Raw%20genetic%20interaction%20datasets:' \
                   '%20Pair-wise%20interaction%20format.zip'
    pmid = 'PMID:27708008'

    # epsilon is an interaction score, see PMID:27708008, supporting material p4
    # "General description of the SGA genetic interaction score"
    effect_type = "epsilon"
    perturbation = "SGA"

    sli_dict = defaultdict(SyntheticLethalInteraction)

    # retrieve Data File S1 and extract and zip SGA_NxN.txt, if necessary
    if not os.path.exists(local_data_file) or force_download:
        with tempfile.TemporaryFile() as temp:
            logging.info("Downloading zip file")
            get_file(temp, zip_file_url)  # download
            logging.info("Extracting interaction data into " + local_data_file)
            with zipfile.ZipFile(temp) as in_zip:
                with in_zip.open(interaction_data_file) as interaction_data,\
                        gzip.open(local_data_file, 'wb') as out_zip:
                    for line in interaction_data.readlines():
                        out_zip.write(line)

    with gzip.open(local_data_file, 'rt') as interaction_data:
        reader = csv.reader(interaction_data, delimiter='\t')
        header = next(reader)
        for row in reader:
            try:
                pvalue = float(row[6])
            except ValueError:
                logging.WARNING("can't convert {} to float".format(row[6]))
            if pvalue > pvalue_cutoff:
                continue

            gene_A_id = "n/a"
            if row[1].upper() in symbol2entrezID:
                gene_A_id = "NCBIGene:{}".format(symbol2entrezID.get(row[1].upper()))

            gene_B_id = "n/a"
            if row[3].upper() in symbol2entrezID:
                gene_B_id = "NCBIGene:{}".format(symbol2entrezID.get(row[3].upper()))

            sli = SyntheticLethalInteraction(gene_A_symbol=row[1],
                                             gene_A_id=gene_A_id,
                                             gene_B_symbol=row[3],
                                             gene_B_id=gene_B_id,
                                             gene_A_pert=perturbation,
                                             gene_B_pert=perturbation,
                                             effect_type=effect_type,
                                             effect_size=row[5],
                                             assay=row[4],
                                             pmid=pmid)

    return defaultdict()


def get_file(file_handle, url):
    """
    Retrieve remote file given fh and url

    :param file_handle: file handle to download file
    :param url: remote url
    """
    with requests.get(url, stream=True) as r:
        total_length = int(r.headers.get('content-length'))
        logging.info("Downloading file from " + url)
        for chunk in progress.bar(r.iter_content(chunk_size=1024),
                                  expected_size=(total_length / 1024) + 1):
            if chunk:
                file_handle.write(chunk)
                file_handle.flush()


symbol2entrezID = EntrezLookup().reverse_lookup
sli_list = parse_luo2009_supplemental_file_S3('data/luo2009.tsv', symbol2entrezID)

for sli in sli_list:
    print(sli.get_tsv_line())

scSymbol2entrezID = EntrezLookup(filename=
                                 "lookup/Saccharomyces_cerevisiae.gene_info.gz",
                                 species_id=["4932", "559292"]
                                 ).reverse_lookup

boone_sli_list = parse_costanzo_boone_2016_NxN_data(symbol2id=scSymbol2entrezID)
for sli in boone_sli_list:
    print(sli.get_tsv_line())

