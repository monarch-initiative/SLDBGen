from enum import Enum


class SlConstants:
    ACTIVATING_MUTATION = 'activating_mutation'
    CELL_VIABILITY_ASSAY = 'cell viability assay'
    CISPLATIN_CYTOTOXICITY_ASSAY = "cisplatin toxicity assay"
    COHORT_STUDY = "cohort study"
    COMPETITIVE_HYBRIDIZATION = 'competitive hybridization'
    CRISPR_CAS9 = "CRISPR CAS9"
    CRISPR_CAS9_INTERFERENCE_ASSAY = "CRISPR-Cas9 Interference assay"
    CYTOTOXICITY_ASSAY = "cytotoxicity assay"
    GROWTH_INHIBITION_ASSAY = "growth inhibition assay"
    KNOCKOUT = 'knockout'
    LOF_MUTATION = 'lof_mutation'
    LOG2_DECREASE_IN_ABUNDANCE = "log2 decrease in abundance"
    MULTICOLOR_COMPETITION_ASSAY = 'multicolor competition assay'
    N_A = "n/a"
    OVEREXPRESSION = "overexpression"
    PHARMACEUTICAL = 'pharmaceutical'
    PHARAMACEUTICAL_INHIBITION_ASSAY = "pharmaceutical inhibition assay"
    RNA_INTERFERENCE_ASSAY = "RNA-interference assay"
    SG_RNA = "sgRNA"
    SH_RNA = 'shRNA'
    SI_RNA = 'siRNA'
    ZSCORE = "Z-score"



    ## Cells

    A549_CELL = "A-549"
    A549_CELLOSAURUS = "CVCL_0023"
    A704_CELL = 'A-704'
    A704_CELLOSAURUS = 'CVCL_1065'
    H358_CELL = 'NCI-H358'
    H358_CELLOSAURUS = 'CVCL_1559'
    HCT_116 ="HCT116"
    HCT_116_CELLOSAURUS = 'CVCL_0291'
    HAP1_CELL = 'HAP1'
    HAP1_CELLOSAURUS = 'CVCL_Y019'
    HELA_CELL = "HeLa-Cells"
    HELA_CELLOSAURUS = "CVCL_0030"
    K562_CELL = "K562 chronic myeloid leukemia cells"
    K562_CELLOSAURUS = "CVCL_0004"
    KBM7_CELL = 'KBM - 7'
    KBM7_CELLOSAURUS = 'CVCL_A426'
    MCF7_CELL = 'MCF-7'
    MCF7_CELLOSAURUS = 'CVCL_0031'
    OVCAR3_CELL = "OVCAR-3"
    OVCAR3_CELLOSAURUS = "CVCL_0465"
    UMRC6_CELL = 'UM-RC-6'
    UMRC6_CELLOSAURUS = 'CVCL_2741'

    ## NCI T
    BREAST_CARCINOMA  = "Breast Carcinoma"
    BREAST_CARCINOMA_NCIT = "NCIT:C4872"
    CHRONIC_MYELOGENOUS_LEUKEMIA = "Chronic Myelogenous Leukemia"
    CHRONIC_MYELOGENOUS_LEUKEMIA_NCIT = "C3174"
    CLEAR_CELL_RENAL_CELL_CARCINOMA = "Clear Cell Renal Cell Carcinoma"
    CLEAR_CELL_RENAL_CELL_CARCINOMA_NCIT = "NCIT:C4033"

    CML_BCRABL_POS = 'Chronic myelogenous leukemia, BCR - ABL1 positive'
    CML_BCRABL_POS_NCIT = 'NCIt: C3174'
    COLON_CARCINOMA = 'Colon carcinoma'
    COLON_CARCINOMA_NCIT = 'NCIT:C4910'
    COLORECTAL_CARCINOMA = "Colorectal Carcinoma"
    COLORECTAL_CARCINOMA_NCIT = "NCIT:C2955"
    GASTRIC_CARCINOMA = 'Gastric Carcinoma'
    GASTRIC_CARCINOMA_NCIT = 'NCIT:C4911'
    LUNG_ADENOCARCINOMA = "Lung adenocarcinoma"
    LUNG_ADENOCARCINOMA_NCIT = "NCIT:C3512"
    MINIMALLY_INVASIVE_LUNG_ADENOCARCINOMA = 'Minimally invasive lung adenocarcinoma'
    MINIMALLY_INVASIVE_LUNG_ADENOCARCINOMA_NCIT = 'NCIT:C2923'
    OVARIAN_SEROUS_ADENOCARCINOMA = "Ovarian serous adenocarcinoma"
    OVARIAN_SEROUS_ADENOCARCINOMA_NCIT ="NCIT:C105555"
    RECURRENT_OVARIAN_CANCER = "Recurrent Ovarian Carcinoma"
    RECURRENT_OVARIAN_CANCER_NCIT = "NCIT:C7833"
    RENAL_CELL_CARCINOMA = ' Renal cell carcinoma'
    RENAL_CELL_CARCINOMA_NCIT = 'NCIT:C9385'




    # Some common gene ids
    ARID1A_GENE_ID = 'NCBIGene:8289'
    BCL2_GENE_ID = 'NCBIGene:596'
    CREBBP_GENE_ID = 'NCBIGene:1387'
    CKS1B_GENE_ID = 'NCBIGene:1163'
    DDR1_GENE_ID = 'NCBIGene:780'
    EP300_GENE_ID = 'NCBIGene:2033'
    EZH2_GENE_ID = 'NCBIGene:2146'
    KRAS_GENE_ID = 'NCBIGene:3845'
    MYC_GENE_ID = 'NCBIGene:4609'
    PBRM1_GENE_ID = 'NCBIGene:55193'
    PDK1_GENE_ID = 'NCBIGene:5163'
    PLK1_GENE_ID = 'NCBIGene:5347'
    PTEN_GENE_ID = 'NCBIGene:5728'
    SLC7A11_GENE_ID = 'NCBIGene:23657'
    SRC_GENE_ID = 'NCBIGene:6714'
    TBK1_GENE_ID = 'NCBIGene:29110'
    TMPRSS4_GENE_ID = 'NCBIGene:56649'
    VHL_GENE_ID = 'NCBIGene:7428'

