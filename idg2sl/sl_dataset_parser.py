import csv
import os
from .parsers.entrez_parser import EntrezParser


class SL_DatasetParser:
    def __init__(self, fname, pmid, symbol2entrezID, hasHeader=True, delim='\t'):
        if not os.path.exists(fname):
            raise ValueError("SL dataset %s does not exist" % fname)
        self.fname = fname
        self.pmid = pmid
        self.symbol2entrezID = symbol2entrezID
        self.rows = []
        with open(fname) as f:
            csvreader = csv.reader(f, delimiter=delim)
            if hasHeader:
                next(csvreader)
            for row in csvreader:
                self.rows.append(row)


    def _mark_maximum_entries(self, sli_dict):
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


    def parse(self):
        """
        This method must be implemented by subclasses
        """
        raise NotImplementedError
