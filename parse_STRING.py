import gzip
import os
from collections import defaultdict





def load_lookup():
    protein2ensembl = defaultdict(str)
    with gzip.open(os.path.join(os.path.dirname(__file__), 'lookup', 'gene2ensembl.gz'), 'rt') as f:
        next(f)  # skip header
        for line in f:
            fields = line.split('\t')
            if fields[0].strip() in ["9606"]:
                ensembl = fields[2]
                if fields[6].startswith("E") and fields[6] is not None:
                    prot_ID = fields[6].split(".")[0]
                    protein2ensembl[prot_ID] = ensembl
    return protein2ensembl


def load_backup():
    protein2ensembl = defaultdict(str)
    with gzip.open(os.path.join(os.path.dirname(__file__), 'lookup', 'prot2gene.gz'), 'rt') as f:
        next(f)  # skip header
        for line in f:
            fields = line.split('\t')
            ensembl = fields[0]
            if fields[1].startswith("E") and fields[1] is not None:
                prot_ID = fields[1]
                protein2ensembl[prot_ID] = ensembl
    return protein2ensembl


def parse_STRING(stringdata, cutoff, outfile, protID2ensembl):
    found = 0
    nfound = 0
    with gzip.open(stringdata, 'rt') as f:
        next(f)  # skip header
        with open(outfile, 'w') as out_f:
            for line in f:
                fields = line.split(' ')
                if int(fields[2]) <= cutoff:
                    continue
                protA = fields[0].split(".")[1]
                protB = fields[1].split(".")[1]
                if protA in protID2ensembl:
                    ensemblA = protID2ensembl.get(protA)
                    found += 1
                #elif protA in protID2ensembl_backup:
                #    ensemblA = protID2ensembl_backup.get(protA)
                #    found += 1
                else:
                    nfound += 1
                    continue
                if protB in protID2ensembl:
                    ensemblB = protID2ensembl.get(protB)
                    found += 1
                #elif protB in protID2ensembl_backup:
                #    ensemblB = protID2ensembl_backup.get(protB)
                #    found += 1
                else:
                    nfound += 1
                    continue
                out_f.write("STRING_data\t" + ensemblA + "\t" + ensemblB + "\t" + fields[2])
    print(f"Found {found} proteins, didn't find {nfound} proteins")
    print(f"{round(nfound / (found + nfound), 2)}% couldn't be found!")


if __name__ == "__main__":
    stringdata = os.path.join(os.path.dirname(__file__), 'data', 'STRING_9606.protein.links.v11.0.txt.gz')
    outfile = "STRING_graph_ensembl.txt"
    cutoff = 700 #threshold on the protein-protein interaction scores

    lookup = load_lookup()      # maps ~93.2% of the IDs
    backup = load_backup()      # maps ~97.9% of the IDs
    # when using both files we also get ~97.9% mapping
    parse_STRING(stringdata, cutoff, outfile, backup)
