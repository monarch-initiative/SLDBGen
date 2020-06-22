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
    assay_string = "RNA-interference assay"
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


def parse_lord_2008(path, symbol2entrezID):
    """
    Parsing data from
    Lord CJ,et al A high-throughput RNA interference screen for DNA repair determinants of PARP inhibitor sensitivity.
    DNA Repair (Amst). 2008 Dec 1;7(12):2010-9. PMID: 18832051.
    This siRNA library targeted over 98% of all of the known DNA repair proteins, as defined by Wood et al. [16] and
    in particular encompassed siRNA targeting all the major components of BER (21 genes), MMR (11 genes),
    NER (28 genes), HR(19 genes) and NHE
    The HTS was performed using the diploid CAL51 breast cancer cell line.
    We used a log 2 surviving fraction (log 2 SF) threshold of −0.1 or less to identify siRNAs that significantly
    sensitized to KU0058948 in the HTS (Supplementary Tables 1–3). This threshold represented the limit of three
    standard deviations from the median of effects in siCON transfected cells. 67 of the 460 experimental siRNA
    satisfied this hit criteria (Fig. 3A).
    Using this distinction, eight genes, plus the control (BRCA1) were identified where both siRNA in the library.
    significantly sensitized to KU0058948 (Fig. 3C); ATR, BRCA2, DDB1, LIG1, PCNA, RAD51, XAB2 and XRCC1.
    We have previously reported, using a different assay system, that silencing of ATR, BRCA2 or RAD51 expression
    significantly sensitizes cells to KU0058948, presumably by causing defective HR [3,6].
    The following parse takes these 8 interactions to be true SL. We additionally generate a list of negative SLs
    using the criteria that neither of the cells attained an effect that was even half as strong
    # The data come as pairs of lines. We will first store the pairs and then only take those that satisfy the
    criteria for positive or negative (PARP sensitization). There are actually 9 SLA pairs including RPA3, which
    is shown in Figure 4C but not mentioned in the text. That is what we get with this parse!
    """
    parp1_symbol = 'PARP1'
    parp1_id = 'NCBIGene:142'
    parp1_perturbation = 'drug'
    gene2_perturbation = 'siRNA'
    pmid = 'PMID:18832051'
    assay_string = "RNA-interference assay"
    effect_type = 'stddev'
    cell_line = 'CAL-51'
    cellosaurus = 'CVCL_1110'
    cancer = "Breast Carcinoma"
    ncit = "NCIT:C4872"
    if not os.path.exists(path):
        raise ValueError("Must path a valid path for Turner et al 2008")
    parpdict = defaultdict(list)
    with open(path) as f:
        # no header
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 3:
                raise ValueError("Malformed line with length %d instead of 3: %s" % (len(fields), line))
            geneBsym = fields[0]
            parp_sens = float(fields[1])
            if geneBsym == 'BRCA1':
                continue
            elif geneBsym == 'GFP-22' or geneBsym == 'SCRAM':
                continue # a control siRNA
            # ignore the third field
            parpdict[geneBsym].append(parp_sens)
    sli_list = []
    for geneBsym, parp_sens_list in parpdict.items():
        if geneBsym in symbol2entrezID:
            geneB_id = "NCBIGene:{}".format(symbol2entrezID.get(geneBsym))
        else:
            geneB_id = "n/a"
        if len(parp_sens_list) != 2:
            raise ValueError("Length of list not equal to 2 for %s (len was %d)" % (geneBsym, len(parp_sens_list))) # should never happen
        if parp_sens_list[0] <= -0.1 and parp_sens_list[1] <= -0.1:
            sli = SyntheticLethalInteraction(gene_A_symbol=parp1_symbol,
                                             gene_A_id=parp1_id,
                                             gene_B_symbol=geneBsym,
                                             gene_B_id=geneB_id,
                                             gene_A_pert=parp1_perturbation,
                                             gene_B_pert=gene2_perturbation,
                                             effect_type=effect_type,
                                             effect_size=min(parp_sens_list[0] ,parp_sens_list[1]),
                                             cell_line=cell_line,
                                             cellosaurus_id=cellosaurus,
                                             cancer_type=cancer,
                                             ncit_id=ncit,
                                             assay=assay_string,
                                             pmid=pmid,
                                             SL=True)
            sli_list.append(sli)
        else:
            effectsize = min(parp_sens_list[0] ,parp_sens_list[1])
            if effectsize > -0.05:
                sli = SyntheticLethalInteraction(gene_A_symbol=parp1_symbol,
                                                 gene_A_id=parp1_id,
                                                 gene_B_symbol=geneBsym,
                                                 gene_B_id=geneB_id,
                                                 gene_A_pert=parp1_perturbation,
                                                 gene_B_pert=gene2_perturbation,
                                                 effect_type=effect_type,
                                                 effect_size=min(parp_sens_list[0], parp_sens_list[1]),
                                                 cell_line=cell_line,
                                                 cellosaurus_id=cellosaurus,
                                                 cancer_type=cancer,
                                                 ncit_id=ncit,
                                                 assay=assay_string,
                                                 pmid=pmid,
                                                 SL=False)
                sli_list.append(sli)
    return sli_list


def parse_toyoshima_2008(path, symbol2entrezID):
    """
    Parsing data from
    Toyoshima M, et al. Functional genomics identifies therapeutic targets for MYC-driven cancer.
    Proc Natl Acad Sci U S A. 2012 Jun 12;109(24):9545-50. PMID: 22623531
    The results of the screen revealed 148 hits, defined according to a Z score of ≥2 (23), including 140 genes and
    eight microRNAs (Fig. 1B). Here, we focus on the 140 gene hits, which we designate MYC-synthetic lethal (MYC-SL)
    genes. To eliminate siRNAs that exhibited substantial growth inhibition properties in normal cells, siRNAs
    with >50% reduced viability in HFF-pBabe were eliminated from further consideration regardless of differential
    toxicity. This process left 102 MYC-SL gene hits for follow-up
    NOTE: I see only 101 genes in the Supplemental table.
    """
    mycsymbol = 'MYC'
    myc_id = 'NCBIGene:4609'
    myc_perturbation = 'overexpression'
    gene2_perturbation = 'siRNA'
    pmid = 'PMID:22623531'
    assay_string = "RNA-interference assay"
    effect_type = 'stddev'
    cell_line = 'HFF-Myc'
    cellosaurus = 'CVCL_Y511'
    cancer = "n/a"
    ncit = "n/a"
    if not os.path.exists(path):
        raise ValueError("Must path a valid path for Turner et al 2008")
    sl_list = []
    with open(path) as f:
        next(f) # skip header
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 6:
                raise ValueError("Bad line with %d fields: %s" % (len(fields), line))
            # Gene Symbol	Accession number	Z score (>than)	%Viability HFF-pB	%Viability HFF-MYC	Ratio pBabe/Myc
            geneBsym = fields[0]
            if geneBsym in symbol2entrezID:
                geneB_id = "NCBIGene:{}".format(symbol2entrezID.get(geneBsym))
            else:
                geneB_id = "n/a"
            accession = fields[1]
            zscore = float(fields[2])
            viability_HFF = float(fields[3])
            viability_HFF_MYC = float(fields[4])
            ratio_pBabe_Myc = float(fields[5])
            sli = SyntheticLethalInteraction(gene_A_symbol=mycsymbol,
                                             gene_A_id=myc_id,
                                             gene_B_symbol=geneBsym,
                                             gene_B_id=geneB_id,
                                             gene_A_pert=myc_perturbation,
                                             gene_B_pert=gene2_perturbation,
                                             effect_type=effect_type,
                                             effect_size=zscore,
                                             cell_line=cell_line,
                                             cellosaurus_id=cellosaurus,
                                             cancer_type=cancer,
                                             ncit_id=ncit,
                                             assay=assay_string,
                                             pmid=pmid,
                                             SL=True)
            sl_list.append(sli)
    return sl_list


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
    ncit = ""			# NCI Thesaurus, Ontology

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
                                     gene_B_symbol = "CSNK2A1",
                                     gene_B_id = "NCBIGene:{}".format(symbol2entrezID.get("CSNK2A1")),
                                     gene_A_pert=gene1_perturbation,
                                     gene_B_pert=gene2_perturbation,
                                     effect_type=effect_type,
                                     effect_size = -0.82,
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
                                     effect_size= -0.96,
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
    cell_line = "Ba/F3"                             # ?
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

            geneB_sym = fields[0]

            if geneB_sym in symbol2entrezID:
                geneB_id = "NCBIGene:{}".format(symbol2entrezID.get(geneB_sym))
            else:
                geneB_id = 'n/a'

            effect = float(fields[8].replace(",", "."))

            threshold = -3                                  # which cutoff?

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

            effect = -4                             # ?

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
                    #effect = float("-2.5")
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
