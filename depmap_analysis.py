from typing import Any, Union
import sys
import wget
import pandas as pd
import os
import numpy as np
from collections import defaultdict
import argparse
from idg2sl.depmap_classes import Gene, SyntheticLethalInteraction
import multiprocessing as mp
import gzip

import time

genes = []

def arg_parser():
    parser = argparse.ArgumentParser(description='Process the DepMap data by downloading the file, '
                                                 'changing HGNC-symbols to EnsemblIDs,'
                                                 'deleting genes with proportion of effective cells above threshold'
                                                 'and blablabla')
    parser.add_argument('--effect', metavar="gene_effect_threshold", type=float, default=-0.5,
                        help='Below which score should a gene be considered effective in a cell line? (default -0.5)')
    parser.add_argument('--prop', metavar="cell_proportion_threshold", type=float, default=0.2,
                        help='Below which proportion of effective cell lines is a gene considered? (default 0.2)')
    parser.add_argument('--intersect', metavar="intersect_threshold", type=float, default=0,
                        help='How many percent of cell lines should be covered by both gene A and gene B? (default 0)')
    args = parser.parse_args()
    return args


def get_data(gene_effect_threshold, cell_proportion_threshold):
    gene_effect_threshold = gene_effect_threshold
    cell_proportion_threshold = cell_proportion_threshold
    global genes

    print("Reducing the original data...")
    url = 'https://ndownloader.figshare.com/files/14221385'
    filename = 'data/gene_effect_corrected.csv'  # path to the file

    if not os.path.exists(filename):
        filename = wget.download(url)

    df = pd.read_csv(filename)  # , index_col=0)
    df = df.drop(df.columns[[0]], axis=1)

    symbol_lookup = defaultdict(str)

    with gzip.open(os.path.join(os.path.dirname(__file__), 'lookup', 'Homo_sapiens.gene_info.gz'), 'rt') as f:
        for line in f:
            fields = line.split('\t')
            if fields[0].strip() in ["9606"]:
                symbol = fields[2]
                ensembl = fields[5].split(":")[-1]
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
            print(f"Data at {i} genes!")
        col = ensembl_df.iloc[:, i]  # corresponds to a column (which shows one gene)
        effective_cells = col.index[col < gene_effect_threshold].tolist()
        proportion = len(effective_cells) / ensembl_df.shape[0]
        i += 1
        if not col.name.startswith("ENS"):
            continue
        if not (0 < proportion < cell_proportion_threshold):
            continue
        gene = Gene(gene_id=col.name, effective_cells=effective_cells)
        genes.append(gene)

    # ensembl_df.to_csv(os.path.join(os.path.dirname(__file__), 'data', 'ensembl_df.csv'), index=False)
    print("Calculation done!")
    # return genes


def intersection(lst1, lst2):
    return [value for value in lst1 if value in set(lst2)]


def f(i, genes, intersect_threshold):
    gene_list = []
    gene_A = genes[i]
    j = i+1
    while j in range(len(genes)):
        gene_B = genes[j]
        cells_covered = len(set(gene_A.get_effective_cells() + gene_B.get_effective_cells()))
        intersect = len(intersection(gene_A.get_effective_cells(), gene_B.get_effective_cells()))
        if intersect / cells_covered > intersect_threshold:
            sli = SyntheticLethalInteraction(gene_A_id=gene_A.get_id(), gene_B_id=gene_B.get_id(), SL=1)
            gene_list.append(sli)
        j += 1
    return gene_list


def get_SLI(intersect_threshold):
    results = []
    sls = []
    p = mp.Pool()
    for i in range(len(genes)):
        results.append(p.starmap_async(f, [(i, genes, intersect_threshold)]))

    #results.append(f(gene.get_id(), gene.get_effective_cells(), genes, intersect_threshold))
    #genes.remove(gene)
    p.close()
    print("Pools are closed!")
    p.join()
    print("Pool joining done")
    print(type(results))
    print(results[0])
    for res in results:
        for list in res.get():
            for sli in list:
                dict = {"gene_A": sli.get_gene_A_id(), "gene_B": sli.get_gene_B_id(), "SL": sli.get_SL()}
                sls.append(dict)
    print("sls erstellt")
    return sls


if __name__ == '__main__':
    start = time.process_time()
    args = arg_parser()

    gene_effect_threshold = args.effect  # threshold from which genes are considered effective
    cell_proportion_threshold = args.prop  # proportion of cells, geneA is considered effective in -> geneB as well!!!
    intersect_threshold = args.intersect  # proportion of cell lines covered by geneA and geneB

    print(f"\ngenes considered effective in cell line below: {args.effect}\n"
          f"necessary proportion of effective cells lines to be considered: {args.prop}\n"
          f"minimum intersect between the set of effective cells from two genes: {args.intersect}\n")

    get_data(gene_effect_threshold, cell_proportion_threshold)
    print("starting the big boys now!")

    sls = get_SLI(intersect_threshold)

    sl_pairs = pd.DataFrame(sls)
    print(sl_pairs.head())
    #sl_pairs.to_csv(f'output_{args.effect}_{args.prop}_{args.intersect}.csv', index=False)

    print(f"\nCalculation took {time.process_time() - start} seconds!\n")
