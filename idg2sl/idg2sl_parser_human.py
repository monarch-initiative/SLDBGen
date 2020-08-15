import logging
import os.path
from _collections import defaultdict
from idg2sl import SyntheticLethalInteraction

from .gene_pair import GenePair

## Some constants
activating_mutation = 'activating_mutation'


def mark_maximum_entries(sli_dict):
    """
    The parsing functions add all SLIs for gene A & B to a list
    Here, we get a dictionary of lists (the list can have one or more entry)
    The keys are GenePair objects.
    We need to mark one entry in each list as being the Max=True
    """
    sli_list = []
    for k, vlist in sli_dict.items():
        vlist.sort(key=lambda x: abs(x.effect_size), reverse=True)
        sli = vlist[0]
        sli.set_maximum()
        sli_list.append(sli)
        for s in vlist[1:]:
            sli_list.append(s)
    return sli_list





def parse_han_2017(path, symbol2entrezID):
    # using supplemental file 14
    gene1_perturbation = "sgRNA"
    gene2_perturbation = "sgRNA"
    pmid = "28319085"
    assay = "RNA Interference assay"
    effect_type = "z-Score"
    cell_line = "K562 chronic myeloid leukemia cells"
    cellosaurus = "CVCL_0004"
    cancer = "Chronic Myelogenous Leukemia"
    ncit = "C3174"

    sli_list = []
    if not os.path.exists(path):
        raise ValueError("Must enter a valid path for Han et al 2017")
    # The following keeps track of the current largest effect size SLI for any given gene A/gene B pair
    sli_dict = defaultdict(list)
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 4:
                logging.error("Only got %d fields but was expecting at least 4 tab-separated fields" % len(fields))

            # seperate genes
            genes = fields[0].split("__")
            geneA_sym = genes[0]
            geneB_sym = genes[1]

            if geneA_sym in symbol2entrezID:
                geneA_id = "NCBIGene:{}".format(symbol2entrezID.get(geneA_sym))
            else:
                geneA_id = 'n/a'

            if geneB_sym in symbol2entrezID:
                geneB_id = "NCBIGene:{}".format(symbol2entrezID.get(geneB_sym))
            else:
                geneB_id = 'n/a'

            effect = -4  # ?

            sli = SyntheticLethalInteraction(gene_A_symbol=geneA_sym,
                                             gene_A_id=geneA_id,
                                             gene_B_symbol=geneB_sym,
                                             gene_B_id=geneB_id,
                                             gene_A_pert=gene1_perturbation,
                                             gene_B_pert=gene2_perturbation,
                                             effect_type=effect_type,
                                             effect_size=effect,
                                             cell_line=cell_line,
                                             cellosaurus_id=cellosaurus,
                                             cancer_type=cancer,
                                             ncit_id=ncit,
                                             assay=assay,
                                             pmid=pmid,
                                             SL=True)
            gene_pair = GenePair(geneA_sym, geneB_sym)
            sli_dict[gene_pair].append(sli)
    sli_list = mark_maximum_entries(sli_dict)
    return sli_list


def parse_shen_2017(path, symbol2entrezID):
    # using supplemental file 3
    # https://static-content.springer.com/esm/art%3A10.1038%2Fnmeth.4225/MediaObjects/41592_2017_BFnmeth4225_MOESM43_ESM.xlsx

    gene1_perturbation = "CRISPR-Cas9 perturbation"
    gene2_perturbation = "CRISPR-Cas9 perturbation"
    pmid = "28319113"
    assay = "CRISPR-Cas9 Interference assay"
    effect_type = "z-Score"
    cancer = ""
    ncit = ""

    sli_list = []
    if not os.path.exists(path):
        raise ValueError("Must enter a valid path for Wang et al 2017")
    # The following keeps track of the current largest effect size SLI for any given gene A/gene B pair
    sli_dict = defaultdict(list)
    with open(path) as f:
        next(f)  # skip header
        next(f)  # skip first row as well
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 7:
                logging.error("Only got %d fields but was expecting at least 7 tab-separated fields" % len(fields))

            geneA_sym = fields[1]
            if geneA_sym in symbol2entrezID:
                geneA_id = "NCBIGene:{}".format(symbol2entrezID.get(geneA_sym))
            else:
                geneA_id = 'n/a'

            geneB_sym = fields[2]
            if geneB_sym in symbol2entrezID:
                geneB_id = "NCBIGene:{}".format(symbol2entrezID.get(geneB_sym))
            else:
                geneB_id = 'n/a'

            if fields[7] == "Synthetic Lethal":
                SL = True
            else:
                SL = False

            cell_line = fields[6].split(",")

            for i in cell_line:
                if i == "293T":
                    cellosaurus = "CVCL_0161"
                    # effect = float("-2.5")
                    effect = float(fields[5].replace(",", "."))
                if i == "HeLa":
                    cellosaurus = "CVCL_0030"
                    effect = float(fields[3].replace(",", "."))
                if i == "A549":
                    cellosaurus = "CVCL_0023"
                    effect = float(fields[4].replace(",", "."))

                sli = SyntheticLethalInteraction(gene_A_symbol=geneA_sym,
                                                 gene_A_id=geneA_id,
                                                 gene_B_symbol=geneB_sym,
                                                 gene_B_id=geneB_id,
                                                 gene_A_pert=gene1_perturbation,
                                                 gene_B_pert=gene2_perturbation,
                                                 effect_type=effect_type,
                                                 effect_size=effect,
                                                 cell_line=i,
                                                 cellosaurus_id=cellosaurus,
                                                 cancer_type=cancer,
                                                 ncit_id=ncit,
                                                 assay=assay,
                                                 pmid=pmid,
                                                 SL=SL)
                gene_pair = GenePair(geneA_sym, geneB_sym)
                sli_dict[gene_pair].append(sli)
    sli_list = mark_maximum_entries(sli_dict)
    return sli_list
