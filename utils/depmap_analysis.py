from typing import Any, Union

import wget
import pandas as pd
import os
import numpy as np
from collections import defaultdict
import concurrent.futures
from pandas import Series, DataFrame
from pandas.io.parsers import TextFileReader
from idg2sl.depmap_classes import Gene

gene_effect_threshold = -0.5  # threshold from which genes are considered effective
cell_proportion_threshold = 0.2  # proportion of cells, geneA is considered effective in -> geneB as well!!!
intersect_threshold = 0.5  # proportion of cell lines covered by geneA and geneB

effective_cells_per_gene = defaultdict(str)
number_effective_cells_per_gene = defaultdict(int)
genes = []  # safes all genes with class entrys

if not os.path.exists(os.path.join(os.path.dirname(__file__), '..', 'data', 'ensembl_df.csv')):
    print("Calculating ensembl_df")
    url = 'https://ndownloader.figshare.com/files/14221385'
    filename = '../data/gene_effect_corrected.csv'  # path to the file

    if not os.path.exists(filename):
        filename = wget.download(url)

    df = pd.read_csv(filename)  # , index_col=0)
    df = df.drop(df.columns[[0]], axis=1)

    symbol_lookup = defaultdict(str)
    with open(os.path.join(os.path.dirname(__file__), '..', 'lookup', 'lookup.txt'), 'r') as f:
        next(f)  # skip header
        for line in f:
            fields = line.split('\t')
            symbol = fields[0]
            ncbi_id = fields[1]
            ensembl = fields[2]
            symbol_lookup[symbol] = ensembl

    new_columns = []
    for i in list(df.columns):
        symbol = i.split(" ")[0]
        ens_ID = symbol_lookup[symbol]
        new_columns.append(ens_ID)

    ensembl_df = df.copy()
    ensembl_df.columns = new_columns

    for i in range(ensembl_df.shape[1]):
        if i == ensembl_df.shape[1]:
            break
        col = ensembl_df.iloc[:, i]  # corresponds to a column (which shows one gene)
        effective_cells = col.loc[col < gene_effect_threshold]
        proportion = len(effective_cells) / ensembl_df.shape[0]
        if not col.name.startswith("ENS"):
            ensembl_df = ensembl_df.drop(col.name, 1)
            continue
        if proportion < cell_proportion_threshold:
            ensembl_df = ensembl_df.drop(col.name, 1)
            continue

    ensembl_df.to_csv("ensembl_df.csv")
else:
    print("ensembl_df found!")
    ensembl_df = pd.read_csv("lookup/ensembl_df.csv", index_col=0)

# get effective cell lines for each gene
for i in range(ensembl_df.shape[1]):
    col = ensembl_df.iloc[:, i]
    effective_cells = col.loc[col < gene_effect_threshold]
    proportion = len(effective_cells) / ensembl_df.shape[0]
    gene = Gene(gene_id=col.name, effective_cells=effective_cells.index, proportion_effective_cells=proportion)
    genes.append(gene)


def intersection(lst1, lst2):
    return [value for value in lst1 if value in set(lst2)]


data = defaultdict(str)


# row_list = []
# build SL dataframe

def f(gene_A):
    row_list = []
    if 0 < gene_A.get_proportion() < cell_proportion_threshold:
        for gene_B in genes:
            if gene_A.get_id() == gene_B.get_id():
                continue
            cells_covered = len(set(gene_A.get_effective_cells() + gene_B.get_effective_cells()))
            intersect = len(intersection(gene_A.get_effective_cells(), gene_B.get_effective_cells()))
            if intersect / cells_covered > intersect_threshold:
                dict = {"gene_A": gene_A.get_id(), "gene_B": gene_B.get_id(), "SL": 1}
                row_list.append(dict)
    return row_list


# print(row_list)

# if __name__ == '__main__':
#	lock = Lock()
#	for num in range(3000):
#	        Process(target = f, args = (lock, num)).start()

sls = []

with concurrent.futures.ProcessPoolExecutor() as executor:
    results = executor.map(f, genes)

    for result in results:
        sls.extend(result)

# processes = []

# for num in range(1000):
#	q = Queue()
#	p = Process(target = f, args = (num, q))
#	p.start()
#	processes.append(p)
#	new_list = q.get()
#	sls.extend(new_list)


# for process in processes:
#	process.join()	

# Pool().close()
# Pool().join()
# print(row_list)


sl_pairs = pd.DataFrame(sls)
print(sl_pairs.head())
sl_pairs.to_csv('output.csv')
