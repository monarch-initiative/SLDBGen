import gzip
import os
from collections import defaultdict

class EntrezParser:
    """
    class to parse the downloaded Homo_sapiens.gene_info.gz file
    """
    def __init__(self, path):
        self.gene_info_path = path
        if not os.path.exists(path):
            raise FileExistsError("Could not find file %s" % path)
        humanSymbol2entrezID = defaultdict(int)
        with gzip.open(self.gene_info_path,'rt') as f:
            for line in f:
                if not line.startswith("9606"):
                    continue # not really us!
                fields = line.rstrip('\n').split('\t')
                geneid = int(fields[1])
                symbol = fields[2]
                humanSymbol2entrezID[symbol] = geneid
        self.symbol2id = humanSymbol2entrezID

    def get_mapping(self):
        return self.symbol2id