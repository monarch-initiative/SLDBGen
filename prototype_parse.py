import logging
import os.path
import sys

import wget
import gzip
from collections import defaultdict

from utils.entrez_lookup import EntrezLookup


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


def parse_costanzo_boone_2016_data() -> defaultdict:
    """
    Costanzo et al. A global genetic interaction network maps a wiring diagram of
    cellular function. Science. 23 Sep 2016: Vol. 353. Issue 6306.
    DOI: 10.1126/science.aaf1420

    This paper describes the results of massive S. cerevisiae SGA experiment.
    The method parses data describing the effect of knocking out all pairwise
    combinations of non-essential genes (SGA_NxN.txt).

    `Description of Supplemental Data is here
    <http://boonelab.ccbr.utoronto.ca/supplement/costanzo2016/>`_

    `This method parses a data file (SGA_NxN.txt) extracted from this zip file:
    <http://boonelab.ccbr.utoronto.ca/supplement/costanzo2016/data_files/
    Data%20File%20S1_Raw%20genetic%20interaction%20datasets:%20Pair-wise%20
    interaction%20format.zip>`_

    :return: defaultdict with SL interactions
    """
    data_file = 'SGA_NxN.txt'
    zip_file = 'http://boonelab.ccbr.utoronto.ca/supplement/costanzo2016/data_files/' \
               'Data%20File%20S1_Raw%20genetic%20interaction%20datasets:' \
               '%20Pair-wise%20interaction%20format.zip'
    return defaultdict()


symbol2entrezID = EntrezLookup().reverse_lookup
sli_list = parse_luo2009_supplemental_file_S3('data/luo2009.tsv', symbol2entrezID)

for sli in sli_list:
    print(sli.get_tsv_line())


