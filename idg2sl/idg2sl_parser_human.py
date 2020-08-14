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




def parse_Shen2015(path, symbol2entrezID):
    # GeneA is always CHEK1 (pharmaceutically inhibited by AZD7762)
    # GeneB is in data/Shen_2015.tsv
    geneA_symbol = 'CHEK1'
    geneA_id = 'NCBIGene:1111'
    geneA_perturbation = 'pharmaceutical'
    gene2_perturbation = 'siRNA'
    pmid = 'PMID:26437225'
    assay = "RNA-interference assay"
    effect_type = "z-Score"
    cell_line = "HeLa-Cells"
    cellosaurus = "CVCL_0030"
    cancer = ""
    ncit = ""  # NCI Thesaurus, Ontology

    sli_list = []
    if not os.path.exists(path):
        raise ValueError("Must path a valid path for Shen et al 2015")
    # The following keeps track of the current largest effect size SLI for any given gene A/gene B pair
    sli_dict = defaultdict(list)
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            fields = line.rstrip('\n').split('\t')
            fields[0] = fields[0].replace(",", ".")
            if len(fields) < 3:
                logging.error("Only got %d fields but was expecting at least 3" % len(fields))
                i = 0
                for fd in fields:
                    print("%d) %s" % (i, fd))
                    raise ValueError("Malformed line, must have at least 3 tab-separated fields")
            geneB_sym = fields[1]
            if geneB_sym in symbol2entrezID:
                geneB_id = "NCBIGene:{}".format(symbol2entrezID.get(geneB_sym))
            else:
                geneB_id = 'n/a'

            effect = float(fields[0])
            sl_genes = ["FZR1", "RAD17", "RFC1", "BLM", "CDC73", "CDC6", "WEE1"]
            if geneB_sym in sl_genes:
                SL = True
            else:
                SL = False
            sli = SyntheticLethalInteraction(gene_A_symbol=geneA_symbol,
                                             gene_A_id=geneA_id,
                                             gene_B_symbol=geneB_sym,
                                             gene_B_id=geneB_id,
                                             gene_A_pert=geneA_perturbation,
                                             gene_B_pert=gene2_perturbation,
                                             effect_type=effect_type,
                                             effect_size=effect,
                                             cell_line=cell_line,
                                             cellosaurus_id=cellosaurus,
                                             cancer_type=cancer,
                                             ncit_id=ncit,
                                             assay=assay,
                                             pmid=pmid,
                                             SL=SL)
            gene_pair = GenePair(geneA_symbol, geneB_sym)
            sli_dict[gene_pair].append(sli)
    sli_list = mark_maximum_entries(sli_dict)
    return sli_list


def parse_pathak_2015(symbol2entrezID):
    # just 2 SL Interactions, hardcoded
    # SRC Gene is proto-oncogene blocked by Dasatinib
    # trying to maximise the Dasatinib sensitivity by SL interaction

    gene1_symbol = 'SRC'
    gene1_id = 'NCBIGene:6714'
    gene1_perturbation = 'pharmaceutical'
    gene2_perturbation = 'cohort study'
    pmid = 'PMID:26437225'
    assay = "pharmaceutical inhibition study"
    effect_type = "correlation"
    cell_line = "n/a"
    cellosaurus = "n/a"
    cancer = "Recurrent Ovarian Carcinoma"
    ncit = "NCIT:C7833"

    sli_list = []
    sli_dict = defaultdict(list)

    sli = SyntheticLethalInteraction(gene_A_symbol=gene1_symbol,
                                     gene_A_id=gene1_id,
                                     gene_B_symbol="CSNK2A1",
                                     gene_B_id="NCBIGene:{}".format(symbol2entrezID.get("CSNK2A1")),
                                     gene_A_pert=gene1_perturbation,
                                     gene_B_pert=gene2_perturbation,
                                     effect_type=effect_type,
                                     effect_size=-0.82,
                                     cell_line=cell_line,
                                     cellosaurus_id=cellosaurus,
                                     cancer_type=cancer,
                                     ncit_id=ncit,
                                     assay=assay,
                                     pmid=pmid,
                                     SL=True)
    gene_pair = GenePair(gene1_symbol, "CSNK2A1")
    sli_dict[gene_pair].append(sli)
    sli = SyntheticLethalInteraction(gene_A_symbol=gene1_symbol,
                                     gene_A_id=gene1_id,
                                     gene_B_symbol="PRKCE",
                                     gene_B_id="NCBIGene:{}".format(symbol2entrezID.get("PRKCE")),
                                     gene_A_pert=gene1_perturbation,
                                     gene_B_pert=gene2_perturbation,
                                     effect_type=effect_type,
                                     effect_size=-0.96,
                                     cell_line=cell_line,
                                     cellosaurus_id=cellosaurus,
                                     cancer_type=cancer,
                                     ncit_id=ncit,
                                     assay=assay,
                                     pmid=pmid,
                                     SL=True)
    gene_pair = GenePair(gene1_symbol, "PRKCE")
    sli_dict[gene_pair].append(sli)
    sli_list = mark_maximum_entries(sli_dict)
    return sli_list


def parse_srivas_2016(path, symbol2entrezID):
    # using the human SL interactions (supplemental file 4 page 2)
    # https://www.cell.com/molecular-cell/fulltext/S1097-2765(16)30280-5?innerTabgraphical_S1097276516302805=#secsectitle0105
    gene1_perturbation = 'pharmaceutical'
    gene2_perturbation = 'natural (is a TSG)'
    pmid = 'PMID:27453043'
    assay = "pharmaceutical + siRNA"
    effect_type = "z-Score"
    cell_line = "HeLa-Cells"
    cellosaurus = "CVCL_0030"
    cancer = ""
    ncit = ""

    sli_list = []
    if not os.path.exists(path):
        raise ValueError("Must enter a valid path for Srivas et al 2016")
    # The following keeps track of the current largest effect size SLI for any given gene A/gene B pair
    sli_dict = defaultdict(list)
    with open(path) as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 4:
                raise ValueError("Only got %d fields but was expecting at least 4 tab-separated fields" % len(fields))

            # seperate col containing multiple genes
            geneA_sym = fields[1].split(",")
            geneB_sym = fields[2]
            if geneB_sym in symbol2entrezID:
                geneB_id = "NCBIGene:{}".format(symbol2entrezID.get(geneB_sym))
            else:
                geneB_id = 'n/a'

            effect = float(fields[3].replace(",", "."))

            for i in geneA_sym:
                if i in symbol2entrezID:
                    geneA_id = "NCBIGene:{}".format(symbol2entrezID.get(i))
                else:
                    geneA_id = 'n/a'

                sli = SyntheticLethalInteraction(gene_A_symbol=i,
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
                gene_pair = GenePair(i, geneB_sym)
                sli_dict[gene_pair].append(sli)
    sli_list = mark_maximum_entries(sli_dict)
    return sli_list


def parse_wang_2017(path, symbol2entrezID):
    # using supplemental file 4
    # https://www.cell.com/cell/fulltext/S0092-8674(17)30061-2?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867417300612%3Fshowall%3Dtrue#secsectitle0275
    gene1_sym = "NRAS"
    gene1_id = "NCBIGene: 4893"
    gene1_perturbation = "mutation"
    gene2_perturbation = "sgRNA"
    pmid = "PMID:28162770"
    assay = "CRISPR-Cas9 Interference assay"
    effect_type = "log2FoldChange"
    cell_line = "Ba/F3"  # ?
    cellosaurus = "CVCL_0161"
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

            geneB_sym = fields[0].upper()

            if geneB_sym in symbol2entrezID:
                geneB_id = "NCBIGene:{}".format(symbol2entrezID.get(geneB_sym))
            else:
                geneB_id = 'n/a'

            effect = float(fields[8].replace(",", "."))
            # if effect < -2.5:
            #    break

            threshold = -3  # which cutoff?

            if effect < threshold:
                SL = True
            else:
                SL = False

            sli = SyntheticLethalInteraction(gene_A_symbol=gene1_sym,
                                             species_id="10090",
                                             gene_A_id=gene1_id,
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
                                             SL=SL)
            gene_pair = GenePair(gene1_sym, geneB_sym)
            sli_dict[gene_pair].append(sli)
    sli_list = mark_maximum_entries(sli_dict)
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
