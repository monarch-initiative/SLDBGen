from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants


class Chakraborty2017Parser(SL_DatasetParser):
    """
    A focused short hairpin RNA library, as well as CRISPR (clustered regularly interspaced short palindromic repeats)/Cas9
    (CRISPR-associated protein 9) and a pharmacological inhibitor => pVHL-defective ccRCC cells are
    hyperdependent on the H3K27 methyltransferase EZH1 for survival.
    Normalized reads from the ‘Late’ (endpoint) samples were then normalized to reads at t=0 to measure changes in the
    shRNA representation in each replicate. The normalized samples were analyzed by  RIGER [29]to compare the relative
    enrichment or depletion of the shRNAs in cells with pVHL relative to cells without pVH.
    Our screens identified many candidate VHL synthetic lethal interactors, including the H3K27 methyltransferase EZH1.
    We prioritized EZH1 because we found that ccRCCs display a distinctive H3K27 methylation profile and because the KDM6A
    (UTX) H3K27 demethylase is a known ccRCC tumor suppressor
    Here, we take all genes with p<0.1 as the authors do.
    We take genes with p>0.4 and p<0.6 as negatives (we do not take those with higher p values because they are enriched
    in the VHL-/- cells)
    We used a script to extract all of the p values for all of the genes and put them in a dictionary (see below)
    """
    def __init__(self, fname=None):
        pmid = "28701475"
        super().__init__(fname=fname, pmid=pmid)

    def get_sli(self, geneB_sym, geneB_id, pval, slstatus):
        # Gene A should be eighter NRAS or KRAS.
        # These genes had activating mutations in the cell lines
        vhl = 'VHL'
        vhlID = SlConstants.VHL_GENE_ID

        sli = SyntheticLethalInteraction(gene_A_symbol=vhl,
                                         species_id="10090",
                                         gene_A_id=vhlID,
                                         gene_B_symbol=geneB_sym,
                                         gene_B_id=geneB_id,
                                         gene_A_pert=SlConstants.LOF_MUTATION,
                                         gene_B_pert=SlConstants.SH_RNA,
                                         effect_type=SlConstants.PVAL,
                                         effect_size=pval,
                                         cell_line=SlConstants.A498_CELL,
                                         cellosaurus_id=SlConstants.A498_CELLOSAURUS,
                                         cancer_type=SlConstants.CLEAR_CELL_RENAL_CELL_CARCINOMA,
                                         ncit_id=SlConstants.CLEAR_CELL_RENAL_CELL_CARCINOMA_NCIT,
                                         assay=SlConstants.SG_RNA_DEPLETION_ASSAY,
                                         pmid=self.pmid,
                                         SL=slstatus)
        return sli

    def parse(self):
        pvals = {"TET3": 0.0042, "KDM4A": 0.0106, "TDO2": 0.0233, "FTO": 0.0242, "P4HTM": 0.0286, "P4HA1": 0.0428,
                 "ADO": 0.0428, "MLL2": 0.0475, "P4HB": 0.051, "IDO1": 0.0511, "KDM6B": 0.0617, "PTGS2": 0.0739,
                 "ASPH": 0.0915, "EZH1": 0.0995, "LEPREL2": 0.1026, "COQ3": 0.1029, "SMYD4": 0.108, "GLOD5": 0.1333,
                 "JMJD4": 0.1336, "PTGS1": 0.1348, "SMYD3": 0.1389, "HR": 0.1389, "MAP3K4": 0.1419, "EHMT1": 0.1535,
                 "RFP": 0.1591, "PRDM9": 0.1657, "TANC2": 0.1672, "HSPBAP1": 0.1759, "KDM4D": 0.1807, "JMJD1C": 0.1912,
                 "FLT1": 0.2013, "DOT1L": 0.2042, "KDM5D": 0.209, "TANC1": 0.2164, "PRMT3": 0.2192, "GLO1": 0.2318,
                 "FBXL19": 0.247, "MCEE": 0.247, "CDO1": 0.2495, "ALKBH1": 0.2495, "PLOD2": 0.27, "KDM4B": 0.27,
                 "KDM5C": 0.2718, "ASPHD1": 0.2774, "SETDB2": 0.2888, "PHYH": 0.2945, "MLL5": 0.3119, "KDM2B": 0.3225,
                 "SMYD1": 0.3239, "PADI4": 0.3377, "PRMT6": 0.353, "PRMT5": 0.357, "TMLHE": 0.3616,
                 "JMJD7-PLA2G4B": 0.3654, "ETHE1": 0.368, "P4HA3": 0.3717, "CDX2": 0.3725, "KDM3B": 0.3761,
                 "PRDM2": 0.3774, "JMJD6": 0.3803, "MINA": 0.3819, "CARM1": 0.3832, "ALKBH3": 0.3924, "ASH2L": 0.3977,
                 "KDM1A": 0.4007, "JMJD5": 0.4285, "PLOD3": 0.4348, "MEIS1": 0.4376, "BCMO1": 0.454, "ALKBH2": 0.4553,
                 "RPE65": 0.4702, "C17orf101": 0.4731, "JARID2": 0.4801, "EZH2": 0.483, "PADI1": 0.4917,
                 "AS3MT": 0.4975, "SETD8": 0.5088, "LUCIFERASE": 0.5137, "C14orf169": 0.5168, "SUV39H2": 0.5182, "KDM4C": 0.5298,
                 "LEPRE1": 0.5339, "C20orf7": 0.5495, "PRMT8": 0.5568, "JMJD8": 0.5575, "BCO2": 0.5732, "BBOX1": 0.5863,
                 "PADI3": 0.5868, "LEPREL1": 0.5912, "ADI1": 0.5996, "SETD7": 0.6035, "KDM5A": 0.6058, "PHYHD1": 0.6113,
                 "GFP": 0.614, "SETD3": 0.6145, "SMYD2": 0.6205, "MLL": 0.6238, "ALKBH6": 0.6449, "UTY": 0.6467,
                 "SUV420H1": 0.647, "ALKBH5": 0.6534, "THUMPD2": 0.6575, "SETD4": 0.6587, "HPDL": 0.662,
                 "EGLN2": 0.6723, "KDM1B": 0.6821, "SHMT1": 0.6862, "SETD1A": 0.7007, "JHDM1D": 0.7107, "EGLN1": 0.7118, "HGD": 0.7139,
                 "SETD2": 0.7211, "ASPHD2": 0.7233, "MLL3": 0.7314, "PRMT2": 0.7334, "HAAO": 0.7353, "PRMT1": 0.7474,
                 "ASH1L": 0.7485, "KDM3A": 0.7502, "MEMO1": 0.7579, "TET1": 0.7673, "TET2": 0.7762, "FOXO3": 0.7809,
                 "SUV39H1": 0.7843, "HPD": 0.7919, "PHF1": 0.7925, "HIF1AN": 0.7996, "ASXL1": 0.8098, "C9orf64": 0.811,
                 "KDM5B": 0.8118, "lacZ": 0.8285, "NSD1": 0.8416, "METTL13": 0.8417, "PADI2": 0.8541, "KDM2A": 0.8648,
                 "OGFOD2": 0.8676, "P4HA2": 0.8752, "PHF2": 0.879, "HNRNPA2B1": 0.8866, "SETDB1": 0.8886,
                 "SETD6": 0.8946, "PLOD1": 0.9137, "OGFOD1": 0.9181, "SUV420H2": 0.9195, "KDM6A": 0.9225,
                 "IDO2": 0.9347, "GLOD4": 0.9464, "EGLN3": 0.9474, "SETD5": 0.9573, "IGBP1": 0.9579, "PHF8": 0.9743, "EHMT2": 0.9855,
                 "PRMT7": 0.9924}
        sli_list = []
        # MLL2 was a symbol for two different genes, I cannot disambiguate
        unclear_gene_symbols = {'MLL2'}
        for geneB, pval in pvals.items():
            if pval < 0.1:
                if geneB in self.entrez_dict:
                    geneBid = self.get_ncbigene_curie(geneB)
                elif geneB in unclear_gene_symbols:
                    continue
                else:
                    raise ValueError("Could not find gene id for %s in Chakraborty 2017" % geneB)
                sli = self.get_sli(geneB_sym=geneB, geneB_id=geneBid, pval=pval, slstatus=True)
                sli_list.append(sli)
            elif 0.4 < pval < 0.6:
                if geneB in self.entrez_dict:
                    geneBid = self.get_ncbigene_curie(geneB)
                else:
                    continue # don't worry, we are looking for negative genes
                sli = self.get_sli(geneB_sym=geneB, geneB_id=geneBid, pval=pval, slstatus=False)
                sli_list.append(sli)
        return sli_list
