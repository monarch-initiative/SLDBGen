from enum import Enum




class SlConstants(Enum):
    ACTIVATING_MUTATION = 1 #'activating_mutation'
    KNOCKOUT = 2
    SI_RNA = 3
    DRUG = 4


    def to_string(self):
        if self.value == 1:
            return 'activating_mutation'
        elif self.value == 2:
            return 'knockout'
        else:
            return "TODO (SlConstants in sl_dataset_parser"