import logging
import os.path
from _collections import defaultdict
from idg2sl import SyntheticLethalInteraction




## Some constants
activating_mutation = 'activating_mutation'


class GenePair:
    """
    We need a class to act as a key in dictionaries with the key being made of gene A and gene B
    because in some of our datasets, both gene A and gene B can vary
    """
    def __init__(self, geneA, geneB):
        self.gene_A = geneA
        self.gene_B = geneB

    def __hash__(self):
        return hash((self.gene_A, self.gene_B))

    def __eq__(self, other):
        return (self.gene_A, self.gene_B) == (other.gene_A, other.gene_B)

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        return not (self == other)



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
    kras_perturbation = activating_mutation
    gene2_perturbation = 'shRNA'
    pmid = 'PMID:19490893'
    assays = ['competitive hybridization', 'multicolor competition assay']
    assay_string = ";".join(assays)
    effect_type = 'stddev'
    cell_line = "DLD-1"
    cellosaurus = "CVCL_0248"
    cancer = "Colorectal Carcinoma"
    ncit = "NCIT:C2955"

    if not os.path.exists(path):
        raise ValueError("Must path a valid path for Luo et al 2009")
    # The following keeps track of the current largest effect size SLI for any given gene A/gene B pair
    sli_dict = defaultdict(list)
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 8:
                logging.error("Only got %d fields but was expecting 8" % len(fields))
                i = 0
                for fd in fields:
                    print("%d) %s" % (i, fd))
                raise ValueError("Malformed line, must have 8 tab-separated fields")
            geneB_sym = fields[0]
            geneB_refSeq = fields[1]
            if geneB_sym in symbol2entrezID:
                geneB_id = "NCBIGene:{}".format(symbol2entrezID.get(geneB_sym))
            else:
                geneB_id = 'n/a'

            stddev = float(fields[5])
            SL = True  # All data in this set is True # TODO CHECK
            sli = SyntheticLethalInteraction(gene_A_symbol=kras_symbol,
                                             gene_A_id=kras_id,
                                             gene_B_symbol=geneB_sym,
                                             gene_B_id=geneB_id,
                                             gene_A_pert=kras_perturbation,
                                             gene_B_pert=gene2_perturbation,
                                             effect_type=effect_type,
                                             effect_size=stddev,
                                             cell_line=cell_line,
                                             cellosaurus_id=cellosaurus,
                                             cancer_type=cancer,
                                             ncit_id=ncit,
                                             assay=assay_string,
                                             pmid=pmid,
                                             SL=SL)
            gene_pair = GenePair(kras_symbol, geneB_sym)
            sli_dict[gene_pair].append(sli)
    sli_list = mark_maximum_entries(sli_dict)
    return sli_list


def parse_bommi_reddi_2008(path, symbol2entrezID):
    """
     Bommi-Reddy A, et al  Kinase requirements in human cells: III. Altered kinase
    requirements in VHL-/- cancer cells detected in a pilot synthetic lethal screen.
    Proc Natl Acad Sci U S A. 2008 Oct 28;105(43):16484-9. PubMed PMID: 18948595.
    The paper describes the use of two cell types: 786-O (Table S4) and RCC4 (Table S5).
    Each of these cells has a loss of function mutation in the tumor suppressor gene VHL.
    The defects in the other gene were induced with a lentiviral vector that produced a shRNA. The
    final results are in Tables S4 and S5.
    The authors require a differential viability above 35% for the 786-O cells or a differential
    viability above 20% for the RCC4 cells. We have copied that values from the supplemental tables
    (which were PDF files) into a TSV file called data/bommi-reddy-2008.tsv. For the output file,
    we demand that at least one shRNA was above threshold, and just report the best outcome per gene.
    The authors showed experimentally that two genes, IRR and HER4, actually are not SL, so we remove these
    by hand
    :return:
    """
    # because of the experiment, geneA is always KRAS. GeneB is in Table S3
    vhl_symbol = 'VHL'
    vhl_id = 'NCBIGene:7428'
    vhl_perturbation = 'lof_mutation'
    gene2_perturbation = 'shRNA'
    pmid = 'PMID:18948595'
    assays = ['competitive hybridization', 'multicolor competition assay']
    assay_string = ";".join(assays)
    effect_type = 'differential_viability'
    cell_786O = "786-0"
    cellosaurus_786O = "CVCL_1051"
    cell_RCC4 = "RCC4"
    cellosaurus_RCC4 = "CVCL_0498"
    cancer = "Clear Cell Renal Cell Carcinoma"
    ncit = "NCIT:C4033"
    if not os.path.exists(path):
        raise ValueError("Must path a valid path for Bommi-Reddy A, et al 2008")
    # The following keeps track of the current largest effect size SLI for any given gene A/gene B pair
    sli_dict = defaultdict(list)
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 4:
                raise ValueError("Only got %d fields but was expecting 4" % len(fields))
            geneB_sym = fields[0]
            if geneB_sym is "IRR" or geneB_sym is "HER4":
                continue
            if geneB_sym in symbol2entrezID:
                geneB_id = "NCBIGene:{}".format(symbol2entrezID.get(geneB_sym))
            else:
                geneB_id = "n/a"
            effect = float(fields[1])
            cell = fields[2]
            if cell == 'RCC4':
                cell_line = cell_RCC4
                cellosaurus = cellosaurus_RCC4
            elif cell == '786-0':
                cell_line = cell_786O
                cellosaurus = cellosaurus_786O
            else:
                raise ValueError("Did not recognize cell type '%s'" % cell)
            table = fields[3]
            assay_string = "differential viability assay {}({})".format(cell, table)
            SL = True  # All data in this set is True # TODO CHECK
            sli = SyntheticLethalInteraction(gene_A_symbol=vhl_symbol,
                                             gene_A_id=vhl_id,
                                             gene_B_symbol=geneB_sym,
                                             gene_B_id=geneB_id,
                                             gene_A_pert=vhl_perturbation,
                                             gene_B_pert=gene2_perturbation,
                                             effect_type=effect_type,
                                             effect_size=effect,
                                             cell_line=cell_line,
                                             cellosaurus_id=cellosaurus,
                                             cancer_type=cancer,
                                             ncit_id=ncit,
                                             assay=assay_string,
                                             pmid=pmid,
                                             SL=SL)
            gene_pair = GenePair(vhl_symbol, geneB_sym)
            sli_dict[gene_pair].append(sli)
    sli_list = mark_maximum_entries(sli_dict)
    return sli_list


def parse_turner_2008(path, symbol2entrezID):
    """
    Parse data from  Turner NC, et al., A synthetic lethal siRNA screen identifying genes mediating
    sensitivity to a PARP inhibitor. EMBO J. 2008 May 7;27(9):1368-77. PubMed PMID: 18388863
    PARP1 was inhibited by the PARP inhibitor KU0058948, and a short interfering RNA library targeting
    779 human protein kinase and kinase assocaited genes was applied.
    :param path:
    :return:
    """
    parp1_symbol = 'PARP1'
    parp1_id = 'NCBIGene:142'
    parp1_perturbation = 'drug'
    gene2_perturbation = 'siRNA'
    pmid = 'PMID:18388863'
    assays = ['competitive hybridization', 'multicolor competition assay']
    assay_string = ";".join(assays)
    effect_type = 'stddev'
    cell_line = 'CAL-51'
    cellosaurus = 'CVCL_1110'
    cancer = "Breast Carcinoma"
    ncit = "NCIT:C4872"
    if not os.path.exists(path):
        raise ValueError("Must path a valid path for Turner et al 2008")
    sli_dict = defaultdict(list)
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            if len(line) < 3:
                raise ValueError("Bad line for Turner et al")
            fields = line.rstrip('\n').split('\t')
            geneB_sym = fields[0]
            if geneB_sym in symbol2entrezID:
                geneB_id = "NCBIGene:{}".format(symbol2entrezID.get(geneB_sym))
            else:
                geneB_id = "n/a"
            zscore = float(fields[1])
            percent_inhib = fields[2]
            if zscore <= -3.0:
                SL = True
            else:
                SL = False
            sli = SyntheticLethalInteraction(gene_A_symbol=parp1_symbol,
                                             gene_A_id=parp1_id,
                                             gene_B_symbol=geneB_sym,
                                             gene_B_id=geneB_id,
                                             gene_A_pert=parp1_perturbation,
                                             gene_B_pert=gene2_perturbation,
                                             effect_type=effect_type,
                                             effect_size=zscore,
                                             cell_line=cell_line,
                                             cellosaurus_id=cellosaurus,
                                             cancer_type=cancer,
                                             ncit_id=ncit,
                                             assay=assay_string,
                                             pmid=pmid,
                                             SL=SL)
            gene_pair = GenePair(parp1_symbol, geneB_sym)
            sli_dict[gene_pair].append(sli)
        sli_list = mark_maximum_entries(sli_dict)
        return sli_list


def parse_steckel_2012(path, symbol2entrezID):
    """
    Steckel M, et al. Determination of synthetic lethal interactions in KRAS oncogene-dependent cancer cells reveals novel
    therapeutic targeting strategies. Cell Res. 2012 Aug;22(8):1227-45.  PubMed PMID: 22613949
    HCT-116 human colon cancer cells and the isogenic derivative, HKE-3, in which the activated, but not the
    normal, KRAS allele has been removed by homologous recombination
    a Delta-Z-score cut-off value of 3.3 was selected to generate a primary hit list of 89 genes
    (∼ 1.2% of the total number of genes screened)
    KRAS itself was placed ninth in this ΔZ-score ranking list, scoring very highly in HCT-116 and weakly in HKE-3 cells
    and thereby serving as an important internal control (Supplementary information, Table S1). Of the remaining 88 genes,
    18 with a Z-score in excess of 2.0 in the HKE-3 cell line alone were eliminated from further analysis, following the
    reasoning that high levels of apoptosis resulting from siRNA-mediated silencing of these genes would likely
    constitute an undesirably strong cytotoxic effect in wild-type KRAS cells (Supplementary information, Figure S2B
    and S2C). Thus, our starting list for further validation comprised 70 candidate genes.
    6159 lines in Harnessing Supp.
    Parse strategy
    delta_zscore >= 3.3 reveals 89 candidates (correct)
    if we remove kras and genes with HKE3_zscore > 2, we get to 70 genes (correct)
    This does not match with the Harnessing paper, which reports 80 SL genes. We will stick with the results
    of the initial screen, i.e., 70 genes, and count the rest as negatives.
    """
    kras_symbol = 'KRAS'
    kras_id = 'NCBIGene:3845'
    kras_perturbation = activating_mutation
    gene2_perturbation = 'siRNA'
    pmid = 'PMID:22613949'
    assay_string = "RNA-interference  assay"
    effect_type = 'stddev'
    cell_line = 'HCT-116'
    cellosaurus = 'CVCL_0291'
    cancer = "Colorectal Carcinoma"
    ncit = "NCIT:C2955"
    if not os.path.exists(path):
        raise ValueError("Must path a valid path for Turner et al 2008")
    sli_dict = defaultdict(list)
    # GeneID	Locus.ID	Accession	HCT-116 Z-score	HKE-3 Z-score	D Z-score
    with open(path) as f:
        next(f) # skip header
        for line in f:
            F = line.rstrip('\n').split('\t')
            if len(F) != 6:
                raise ValueError("Line has %d fields (should have 6): %s" % (len(F), line))
            geneB_sym = F[0]
            locusID = F[1]
            accession = F[2]
            HCT116_zscore = float(F[3])
            HKE3_zscore = float(F[4])
            delta_zscore = float(F[5])
            if geneB_sym in symbol2entrezID:
                geneB_id = "NCBIGene:{}".format(symbol2entrezID.get(geneB_sym))
            else:
                geneB_id = "n/a"
            if geneB_sym == "KRAS":
                continue # This was an internal control!
            if delta_zscore >= 3.3 and HKE3_zscore < 2:
                SL = True
            else:
                SL = False
            sli = SyntheticLethalInteraction(gene_A_symbol=kras_symbol,
                                             gene_A_id=kras_id,
                                             gene_B_symbol=geneB_sym,
                                             gene_B_id=geneB_id,
                                             gene_A_pert=kras_perturbation,
                                             gene_B_pert=gene2_perturbation,
                                             effect_type=effect_type,
                                             effect_size=HCT116_zscore,
                                             cell_line=cell_line,
                                             cellosaurus_id=cellosaurus,
                                             cancer_type=cancer,
                                             ncit_id=ncit,
                                             assay=assay_string,
                                             pmid=pmid,
                                             SL=SL)
            gene_pair = GenePair(kras_symbol, geneB_sym)
            sli_dict[gene_pair].append(sli)
        sli_list = mark_maximum_entries(sli_dict)
        return sli_list