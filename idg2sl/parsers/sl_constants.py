from enum import Enum




class SlConstants(Enum):
    ACTIVATING_MUTATION = 1 #'activating_mutation'
    KNOCKOUT = 2
    SI_RNA = 3
    DRUG = 4
    LOF_MUTATION = 5
    SH_RNA = 6
    COMPETITIVE_HYBRIDIZATION = 7
    MULTICOLOR_COMPETITION_ASSAY = 8
    RNA_INTERFERENCE_ASSAY = 9
    OVEREXPRESSION = 10



    def to_string(self):
        if self.value == 1:
            return 'activating_mutation'
        elif self.value == 2:
            return 'knockout'
        elif self.value == 3:
            return 'siRNA'
        elif self.value == 4:
            return 'drug'
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
        else:
            return "TODO (SlConstants in sl_dataset_parser"