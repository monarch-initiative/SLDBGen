from enum import Enum


class SlConstants:
    ACTIVATING_MUTATION = 'activating_mutation'
    CELL_VIABILITY_ASSAY = 'cell viability assay'
    CISPLATIN_CYTOTOXICITY_ASSAY = "cisplatin toxicity assay"
    COHORT_STUDY = "cohort study"
    COMPETITIVE_HYBRIDIZATION = 'competitive hybridization'
    CRISPR_CAS9 = "CRISPR CAS9"
    CRISPR_CAS9_INTERFERENCE_ASSAY = "CRISPR-Cas9 Interference assay"
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
    HCT_116 ="HCT116"
    HCT_116_CELLOSAURUS = 'CVCL_0291'
    HAP1_CELL = 'HAP1'
    HAP1_CELLOSAURUS = 'CVCL_Y019'
    KBM7_CELL = 'KBM - 7'
    KBM7_CELLOSAURUS = 'CVCL_A426'

    ## NCI T
    COLON_CARCINOMA = 'Colon carcinoma'
    COLON_CARCINOMA_NCIT = 'NCIT:C4910'
    GASTRIC_CARCINOMA = 'Gastric Carcinoma'
    GASTRIC_CARCINOMA_NCIT = 'NCIT:C4911'
    LUNG_ADENOCARCINOMA = "Lung adenocarcinoma"
    LUNG_ADENOCARCINOMA_NCIT = "NCIT:C3512"
    CML_BCRABL_POS = 'Chronic myelogenous leukemia, BCR - ABL1 positive'
    CML_BCRABL_POS_NCIT = 'NCIt: C3174'

    # Some common gene ids
    ARID1A_GENE_ID = 'NCBIGene:8289'
    CREBBP_GENE_ID = 'NCBIGene:1387'
    EP300_GENE_ID = 'NCBIGene:2033'
    EZH2_GENE_ID = 'NCBIGene:2146'
    KRAS_GENE_ID = 'NCBIGene:3845'
    PTEN_GENE_ID = 'NCBIGene:5728'

