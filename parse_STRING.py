import gzip
import os
import wget
from collections import defaultdict
from utils.lookup import Lookup
import idg2sl

protID2ensembl = Lookup().protein2ensembl
#print(protID2ensembl)
#protID2ensembl_backup = EnsemblLookup().backup_protein_lookup
protID2ensembl_backup = defaultdict(str)
#print(protID2ensembl_backup)

stringdata = os.path.join(os.path.dirname(__file__), 'data', 'STRING_9606.protein.links.v11.0.txt.gz')
outfile = "STRING_graph_ensembl.txt"


def parse_STRING(stringdata, outfile):
    found = 0
    nfound = 0
    with gzip.open(stringdata, 'rt') as f:
        next(f)  # skip header
        with open(outfile, 'w') as out_f:
            for line in f:
                fields = line.split(' ')
                protA = fields[0].split(".")[1]
                protB = fields[1].split(".")[1]
                ensemblA, ensemblB = "", ""
                if protA in protID2ensembl:
                    ensemblA = protID2ensembl.get(protA)
                    found += 1
                elif protA in protID2ensembl_backup:
                    ensemblA = protID2ensembl.get(protA)
                    found += 1
                else:
                    ensemblA = "n/a"
                    nfound += 1
                if protB in protID2ensembl:
                    ensemblB = protID2ensembl.get(protB)
                    found += 1
                elif protB in protID2ensembl_backup:
                    ensemblB = protID2ensembl.get(protB)
                    found += 1
                else:
                    ensemblB = "n/a"
                    nfound += 1
                out_f.write("STRING_data\t" + ensemblA + "\t" + ensemblB + "\t" + fields[2])
    print("Found %d proteins, didn't find %d proteins" % (found, nfound))
    print(nfound/(found + nfound))


parse_STRING(stringdata, outfile)
