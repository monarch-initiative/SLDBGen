## This script generates a simple summary of the SLI data.
## The script `parse_human_SLI.py` needs to be run first

import os
import csv
from collections import defaultdict
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

sl_data_file = 'SL_data.tsv'
if not os.path.isfile(sl_data_file):
    raise ValueError("Could not find %s. Please run parse_human_SLI.py before this script")

all_symbols = set()
all_symbols_positive_set = set()
connections = defaultdict(set)
degrees = defaultdict(int)
sli_count = 0

with open(sl_data_file) as f:
    creader = csv.DictReader(f, delimiter="\t")
    for row in creader: 
        geneA = row['geneA']
        geneB = row['geneB']
        all_symbols.add(geneA)
        all_symbols.add(geneB)
        SL = row['SL']
        if SL == 'T':
            all_symbols_positive_set.add(geneA)
            all_symbols_positive_set.add(geneB)
            connections[geneA].add(geneB)
            connections[geneB].add(geneA)
            sli_count += 1

print("[INFO] Total genes: %d" % len(all_symbols))
print("[INFO] Total genes in positive set: %d" % len(all_symbols_positive_set))
for k, v in connections.items():
    gene = k
    degree = len(v)
    if degree > 50: 
        print("%s - degree: %d" % (gene, degree))
    degrees[degree] += 1
print("[INFO] Degree distribution:")
for k, v in degrees.items():
    print("Degree: %d, number of genes: %d" % (k, v))
print("[INFO] Total number of SLI %d" % sli_count)

# Plot with seaborn. First put the data into a Pandas DataFrame
dict_list = []
degree_one = 0


X = []
Y = []

total = 0
for k, v in degrees.items():
    total += v

for k, v in degrees.items():
    fraction = v/total
    X.append(k)
    Y.append(fraction)

plt.loglog(X,Y,color='b', marker='o', linestyle='None')
plt.xlabel("degree")
plt.ylabel("fraction of nodes")
plt.title("Degree distribution of SLI network")

plt.savefig('sliDegreeDistribution.pdf')

