from typing import Any, Union

import wget
import pandas as pd
import os
import numpy as np
from collections import defaultdict
import concurrent.futures
from pandas import Series, DataFrame
from idg2sl.depmap_classes import Gene
import multiprocessing as mp



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
    row_list = []
    for gene_B in genes:
        if id == gene_B.get_id():
            continue
        cells_covered = len(set(effective_cells + gene_B.get_effective_cells()))
        intersect = len(intersection(effective_cells, gene_B.get_effective_cells()))
        if intersect / cells_covered > intersect_threshold:
            dict = {"gene_A": id, "gene_B": gene_B.get_id(), "SL": 1}
            row_list.append(dict)
    # queue.put(row_list)
    # sls.extend(row_list)
    return row_list


def start_procs():
    results = []
    procs = []
    # queue = multiprocessing.Queue()

    # instantiating process with arguments
    with mp.Pool(4) as p:
        for gene in genes:
            results.append(p.apply(f, args=(gene.get_id(), gene.get_effective_cells())))
        p.close()
        p.join()

    for result in results:
        sls.append(result.get())

    print("all proc started")
    # while not queue.empty():
    #     results.append(queue.get())



def normal():
    i = 0
    for gene in genes:
        i += 1
        if i % 100 == 0:
            print(f"{i} genes done!")
        result = f(gene.get_id(), gene.get_effective_cells())
        sls.append(result)


def log_result(result):
    # This is called whenever f(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    sls.append(result)

def apply_async_with_callback():
    with mp.Pool() as p:
        i = 0
        for gene in genes:
            if i%100 == 0:
                print(i)
            i += 1
            sls.append(p.apply_async(f, (gene.get_id(), gene.get_effective_cells())).get())
    p.close()
    p.join()
    print(sls)





if __name__ == '__main__':
    get_ensembl_df()
    print("starting the big boys now!")
    # apply_async_with_callback()

    results = []
    p =  mp.Pool()
    i = 0
    for gene in genes:
        if i%100 == 0:
            print(i)
        i += 1
        results.append(p.apply_async(f, (gene.get_id(), gene.get_effective_cells(), genes)))
    p.close()
    print("Processes done")
    p.join()
    print("join done")

    sls = [res.get() for res in results]
    print("sls erstellt")

    print(sls[1])

    sl_pairs = pd.DataFrame(sls)
    print(sl_pairs.head())
    sl_pairs.to_csv('output.csv')




































