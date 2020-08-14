from enum import Enum




class SlConstants(Enum):
    ACTIVATING_MUTATION = 1 #'activating_mutation'


    def to_string(self):
        if self.value == 1:
            return 'activating_mutation'
        else:
            return "TODO (SlConstants in sl_dataset_parser"