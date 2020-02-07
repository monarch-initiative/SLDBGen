import csv
import gzip
import logging
import os.path
import tempfile
import zipfile

import requests
from clint.textui import progress
from idg2sl import SyntheticLethalInteraction


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


def parse_costanzo_boone_2016_NxN_data(symbol2id,
                                       force_download=False,
                                       pvalue_cutoff=0.05
                                       ) -> list:
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

    This method sets the SL = True for SL pairs where pvalue < pvalue_cutoff. Other
    SL pairs are still ingested, SL = False. Downstream user can filter further using
    epsilon score. See methods for details:
        "The interaction data ... should be filtered prior to use. We suggest three
        different thresholds [lenient (P < 0.05), intermediate (P < 0.05 and e < 0.08),
        and stringent confidence (P <0.05 and e > 0.16 or e < -0.12)]"

    `Detailed description of Supplemental Data is here:
    <http://boonelab.ccbr.utoronto.ca/supplement/costanzo2016/>`_

    :param symbol2id: dict produced by `EntrezLookup` to look up IDs from symbols
    :param force_download: force dl of zip file and rewriting interaction data [false]
    :param pvalue_cutoff: only SLIs with pvalues less than this will be output [0.05]
    :return: list with SL interactions
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

    sli_list = list()

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
            is_sl = False
            try:
                pvalue = float(row[6])
                if pvalue < pvalue_cutoff:
                    is_sl = True
            except ValueError:
                logging.WARNING("can't convert {} to float".format(row[6]))

            gene_A_id = "n/a"
            if row[1].upper() in symbol2id:
                gene_A_id = "NCBIGene:{}".format(symbol2id.get(row[1].upper()))

            gene_B_id = "n/a"
            if row[3].upper() in symbol2id:
                gene_B_id = "NCBIGene:{}".format(symbol2id.get(row[3].upper()))

            sli = SyntheticLethalInteraction(gene_A_symbol=row[1],
                                             gene_A_id=gene_A_id,
                                             gene_B_symbol=row[3],
                                             gene_B_id=gene_B_id,
                                             gene_A_pert=perturbation,
                                             gene_B_pert=perturbation,
                                             effect_type=effect_type,
                                             effect_size=row[5],
                                             assay=row[4],
                                             pmid=pmid,
                                             SL=is_sl)
            sli_list.append(sli)
    return sli_list
