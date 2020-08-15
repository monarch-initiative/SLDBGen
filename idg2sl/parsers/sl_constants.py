from enum import Enum




class SlConstants(Enum):
    ACTIVATING_MUTATION = 1 #'activating_mutation'
    KNOCKOUT = 2
    SI_RNA = 3
    PHARMACEUTICAL = 4
    LOF_MUTATION = 5
    SH_RNA = 6
    COMPETITIVE_HYBRIDIZATION = 7
    MULTICOLOR_COMPETITION_ASSAY = 8
    RNA_INTERFERENCE_ASSAY = 9
    OVEREXPRESSION = 10
    COHORT_STUDY = 11
    PHARAMACEUTICAL_INHIBITION_ASSAY = 12
    GROWTH_INHIBITION_ASSAY = 13
    CRISPR_CAS9_INTERFERENCE_ASSAY = 14
    SG_RNA = 15
    CRISPR_CAS9 = 16
    ZSCORE = 17




    def to_string(self):
        if self.value == 1:
            return 'activating_mutation'
        elif self.value == 2:
            return 'knockout'
        elif self.value == 3:
            return 'siRNA'
        elif self.value == 4:
            return 'pharmaceutical'
        elif self.value == 5:
            return 'lof_mutation'
        elif self.value == 6:
            return 'shRNA'
        elif self.value == 7:
            return 'competitive hybridization'
        elif self.value == 8:
            return 'multicolor competition assay'
        elif self.value == 9:
            return "RNA-interference assay"
        elif self.value == 10:
            return "overexpression"
        elif self.value == 11:
            return "cohort study"
        elif self.value == 12:
            return "pharmaceutical inhibition assay"
        elif self.value == 13:
            return "growth inhibition assay"
        elif self.value == 14:
            return "CRISPR-Cas9 Interference assay"
        elif self.value == 15:
            return "sgRNA"
        elif self.value == 16:
            return "CRISPR CAS9"
        elif self.value == 17:
            return "Z-score"
        else:
            return "TODO (SlConstants in sl_dataset_parser"