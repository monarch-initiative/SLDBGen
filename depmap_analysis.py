from typing import Any, Union

import wget
import pandas as pd
import os
import numpy as np
from collections import defaultdict
import concurrent.futures
from idg2sl.depmap_classes import Gene, SyntheticLethalInteraction
import multiprocessing as mp
import time

start = time.process_time()

gene_effect_threshold = -0.5  # threshold from which genes are considered effective
cell_proportion_threshold = 0.2  # proportion of cells, geneA is considered effective in -> geneB as well!!!
intersect_threshold = 0  # proportion of cell lines covered by geneA and geneB

effective_cells_per_gene = defaultdict(str)
number_effective_cells_per_gene = defaultdict(int)
sls = []    # list containing final SL pairs
genes = []  # safes all genes with class entrys


def get_ensembl_df():

    if not os.path.exists(os.path.join(os.path.dirname(__file__), 'data', 'ensembl_df.csv')):
        print("Calculating ensembl_df")
        url = 'https://ndownloader.figshare.com/files/14221385'
        filename = 'data/gene_effect_corrected.csv'  # path to the file

        if not os.path.exists(filename):
            filename = wget.download(url)

        df = pd.read_csv(filename)  # , index_col=0)
        df = df.drop(df.columns[[0]], axis=1)

        symbol_lookup = defaultdict(str)
        with open(os.path.join(os.path.dirname(__file__), 'lookup', 'lookup.txt'), 'r') as f:
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
        i = 0
        while i in range(ensembl_df.shape[1]):
            if i % 1000 == 0:
                print(f"{i} genes analysed")
            col = ensembl_df.iloc[:, i]  # corresponds to a column (which shows one gene)
            effective_cells = col.index[col < gene_effect_threshold].tolist()
            proportion = len(effective_cells) / ensembl_df.shape[0]
            if not col.name.startswith("ENS"):
                ensembl_df = ensembl_df.drop(col.name, 1)
                continue
            if not (0 < proportion < cell_proportion_threshold):
                ensembl_df = ensembl_df.drop(col.name, 1)
                continue
            i += 1
            gene = Gene(gene_id=col.name, effective_cells=effective_cells)
            genes.append(gene)

        print(len(genes))

        ensembl_df.to_csv(os.path.join(os.path.dirname(__file__), 'data', 'ensembl_df.csv'))
        print("Done calculationg ensembl_df!")
        print(ensembl_df.head())
    else:
        print("ensembl_df found!")
        ensembl_df = pd.read_csv(os.path.join(os.path.dirname(__file__), 'data', 'ensembl_df.csv'), index_col=0)
        print(ensembl_df.head())
        # get effective cell lines for each gene
        for i in range(ensembl_df.shape[1]):
            col = ensembl_df.iloc[:, i]  # corresponds to a column (which shows one gene)
            effective_cells = col.index[col < gene_effect_threshold].tolist()
            gene = Gene(gene_id=col.name, effective_cells=effective_cells)
            genes.append(gene)


def intersection(lst1, lst2):
    return [value for value in lst1 if value in set(lst2)]


def f(id, effective_cells, genes):
    gene_list = []
    for gene_B in genes:
        if id == gene_B.get_id():
            continue
        cells_covered = len(set(effective_cells + gene_B.get_effective_cells()))
        intersect = len(intersection(effective_cells, gene_B.get_effective_cells()))
        if intersect / cells_covered > intersect_threshold:
            sli = SyntheticLethalInteraction(gene_A_id=id, gene_B_id=gene_B.get_id(), SL=1)
            #sls.append(sli)
            #dict = {"gene_A": id, "gene_B": gene_B.get_id(), "SL": 1}
            gene_list.append(sli)
    # queue.put(row_list)
    # sls.extend(row_list)
    return gene_list




if __name__ == '__main__':
    get_ensembl_df()
    print("starting the big boys now!")
    # apply_async_with_callback()

    results = []
    p = mp.Pool()
    i = 0
    for gene in genes:
        if i % 100 == 0:
            print(i)
        i += 1
        # results.append(p.apply_async(f, (gene.get_id(), gene.get_effective_cells(), genes)))
        results.append(p.starmap_async(f, [(gene.get_id(), gene.get_effective_cells(), genes)]))
    p.close()
    print("Pool closing done")
    p.join()
    print("Pool joining done")
    for res in results:
        for list in res.get():
            for sli in list:
                dict = {"gene_A": sli.get_gene_A_id(), "gene_B": sli.get_gene_B_id(), "SL": sli.get_SL()}
                sls.append(dict)
    print("sls erstellt")



    sl_pairs = pd.DataFrame(sls)
    print(sl_pairs.head())
    sl_pairs.to_csv('output.csv')
    print(f"\nCalculation took {time.process_time() - start} seconds!")




































