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
    genes = []

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
        if ensembl_df.shape[1] % 1000 == 0:
            print(f"Data shortened to {ensembl_df.shape[1]} genes!")
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

    #ensembl_df.to_csv(os.path.join(os.path.dirname(__file__), 'data', 'ensembl_df.csv'), index=False)
    print("Calculation done!")
    return ensembl_df, genes


def intersection(lst1, lst2):
    return [value for value in lst1 if value in set(lst2)]


def f(id, effective_cells, genes, intersect_threshold):
    gene_list = []
    for gene_B in genes:
        if id == gene_B.get_id():
            continue
        cells_covered = len(set(effective_cells + gene_B.get_effective_cells()))
        intersect = len(intersection(effective_cells, gene_B.get_effective_cells()))
        if intersect / cells_covered > intersect_threshold:
            sli = SyntheticLethalInteraction(gene_A_id=id, gene_B_id=gene_B.get_id(), SL=1)
            gene_list.append(sli)
    return gene_list


if __name__ == '__main__':
    start = time.process_time()
    args = arg_parser()

    gene_effect_threshold = args.effect  # threshold from which genes are considered effective
    cell_proportion_threshold = args.prop  # proportion of cells, geneA is considered effective in -> geneB as well!!!
    intersect_threshold = args.intersect  # proportion of cell lines covered by geneA and geneB

    print(f"\ngenes considered effective in cell line below: {args.effect}\n"
          f"necessary proportion of effective cells lines to be considered: {args.prop}\n"
          f"minimum intersect between the set of effective cells from two genes: {args.intersect}\n")

    effective_cells_per_gene = defaultdict(str)
    number_effective_cells_per_gene = defaultdict(int)
    sls = []  # list containing final SL pairs
    genes = []  # safes all genes with class entries

    ensembl_df, genes = get_data(gene_effect_threshold, cell_proportion_threshold)
    print(ensembl_df.head())
    print("starting the big boys now!")

    results = []
    p = mp.Pool()
    i = 0
    for gene in genes:
        if i % 100 == 0:
            print(i)
        i += 1
        # results.append(p.apply_async(f, (gene.get_id(), gene.get_effective_cells(), genes)))
        results.append(p.starmap_async(f, [(gene.get_id(), gene.get_effective_cells(), genes, intersect_threshold)]))
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
    #sl_pairs.to_csv('output.csv', index=False)

    print(f"\nCalculation took {time.process_time() - start} seconds!\n")




































