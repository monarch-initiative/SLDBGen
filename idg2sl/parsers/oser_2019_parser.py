from idg2sl import SyntheticLethalInteraction
from idg2sl.sl_dataset_parser import SL_DatasetParser
from .sl_constants import SlConstants


class Oser2019Parser(SL_DatasetParser):
    """
    To identify synthetic lethal interactors with RB1 in SCLC, the authors first infected two RB1−/− SCLC cell lines
    (NCI-H82 and NCI-H69) with a lentivirus that expresses RB1 in the presence of doxycycline (DOX-On RB1) or with the
    corresponding empty vector (DOX-On EV). Given that RB1 reexpression had no gross effect on cell proliferation in
    NCI-H82 cells, we used the DOX-On RB1 NCI-H82 cells to perform a RB1 synthetic lethal screen using CRISPR/Cas9-based
    gene editing. Using this strategy, we also identified 104 genes that were synthetic lethal with RB1, based on sgRNA
     depletion in the RB1-deficient cells (DOX-On RB1; no DOX) compared with the RB1-proficient cells (DOX-On RB1; +DOX),
     using a P value cutoff of <0.05 (Fig. 1G; Supplementary Tables S2 and S4).
     The authors focused on AURKB because it was the highest-scoring “druggable” hit.
    """
    def __init__(self, fname=None):
        pmid = '30373918'
        super().__init__(fname=fname, pmid=pmid)

    def parse(self):
        # The following genes were extracted from Supplementary table 4.
        # the value is the p-value
        # 103 genes were obtained with p<0.05 (the manuscript says there are 104)
        oserGenes = {'SNRPG': '0.004200', 'POLR2B': '0.022800', 'POLR2D': '0.000200', 'CCNH': '0.000700',
                     'POLR2A': '0.000900', 'AURKB': '0.001100', 'MRGBP': '0.001500', 'SNRPF': '0.001800',
                     'POLR2G': '0.002500', 'RBM14': '0.003100', 'UHRF1': '0.003500', 'DDB1': '0.004100',
                     'RAN': '0.005000', 'U2AF1': '0.008500', 'POLR2L': '0.018100', 'TADA2A': '0.018800',
                     'UBE2I': '0.021600', 'TTF2': '0.034000', 'TEX10': '0.037600', 'NACA': '0.000100',
                     'PLK1': '0.000100', 'GTF2A2': '0.000700', 'GTF3C1': '0.000900', 'FBL': '0.001100',
                     'ACTL6A': '0.001200', 'TONSL': '0.001200', 'MED8': '0.001400', 'H3F3A': '0.001500',
                     'PELP1': '0.002000', 'SMC4': '0.002500', 'RNMT': '0.002700', 'CDAN1': '0.003200',
                     'HCFC1': '0.003500', 'HELZ2': '0.004800', 'CDK7': '0.004900', 'DDX54': '0.005800',
                     'GTF2B': '0.006700', 'SRSF11': '0.007000', 'SNW1': '0.010000', 'SNAPC4': '0.010100',
                     'CPSF2': '0.010500', 'PAXIP1': '0.010700', 'POLR2C': '0.011200', 'KAT7': '0.011400',
                     'UPF1': '0.012500', 'GTF2H2C': '0.013100', 'C1QBP': '0.013200', 'UBTF': '0.013400',
                     'OGDHL': '0.013600', 'IPO9': '0.013600', 'CLP1': '0.014000', 'DIDO1': '0.014500',
                     'CASP8AP2': '0.015500', 'VPRBP': '0.015600', 'BRF1': '0.016700', 'PPP4R2': '0.016800',
                     'INCENP': '0.017100', 'WDHD1': '0.017200', 'PER2': '0.017900', 'GTF3C3': '0.018000',
                     'NAA50': '0.018100', 'NCAPG': '0.019700', 'U2AF2': '0.019800', 'PDS5A': '0.020700',
                     'SSRP1': '0.021900', 'MTF2': '0.022300', 'NPM1': '0.023700', 'TICRR': '0.023800',
                     'ZNF593': '0.024900', 'HIST1H2BN': '0.025200', 'HMGB1': '0.026400', 'CTCF': '0.026600',
                     'GTF3C5': '0.026800', 'NELFCD': '0.027400', 'MIS18BP1': '0.028300', 'POLR3A': '0.028500',
                     'MED11': '0.028900', 'UBE2N': '0.029700', 'SMC2': '0.030200', 'CTDP1': '0.033300',
                     'ASH2L': '0.033400', 'RNF168': '0.034400', 'WNT4': '0.034500', 'CHD8': '0.035700',
                     'PRDM10': '0.035800', 'PCGF6': '0.036100', 'PRMT1': '0.037300', 'PRPF31': '0.037400',
                     'HMGB3': '0.037800', 'GTF2H4': '0.038400', 'PRDM7': '0.038900', 'AIFM1': '0.040400',
                     'MED24': '0.040500', 'UTP3': '0.041700', 'DR1': '0.041900', 'NCBP1': '0.043100',
                     'PCF11': '0.043200', 'WDR82': '0.045000', 'RBBP4': '0.046100', 'ESPL1': '0.047500',
                     'PMF1': '0.048000', 'CENPP': '0.048800', 'HIST1H2BL': '0.047800'}
        rb1 = 'RB1'
        geneAid = self.get_ncbigene_curie(rb1)
        sli_list = []
        for k, v in oserGenes.items():
            geneB = self.get_current_symbol(k)
            if geneB in self.entrez_dict:
                geneBid = self.get_ncbigene_curie(geneB)
            else:
                raise ValueError("Could not find id for %s in Oser 2019" % geneB)
            sli = SyntheticLethalInteraction(gene_A_symbol=rb1,
                                             gene_A_id=geneAid,
                                             gene_B_symbol=geneB,
                                             gene_B_id=geneBid,
                                             gene_A_pert=SlConstants.LOF_MUTATION,
                                             gene_B_pert=SlConstants.CRISPR_CAS9,
                                             effect_type=SlConstants.PVAL,
                                             effect_size=v,
                                             cell_line=SlConstants.H82_CELL,
                                             cellosaurus_id=SlConstants.H82_CELLOSAURUS,
                                             cancer_type=SlConstants.N_A,
                                             ncit_id=SlConstants.N_A,
                                             assay=SlConstants.SG_RNA_DEPLETION_ASSAY,
                                             pmid=self.pmid,
                                             SL=True)
            sli_list.append(sli)
        return sli_list
