#!/usr/bin/env python3
# coding: utf-8



import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import re
import os
import sys
import argparse
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
sc.settings.verbosity = 3
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=300)


# # Pipeline



def ReadData(feature, phen):
    RawData = ad.read_csv(feature, delimiter='\t', first_column_names=True)
    phen = pd.read_csv(phen, sep='\t', header=0, index_col=0)
    phen.index = [str(i) for i in phen.index]
    RawData.obs['phen'] = phen.loc[RawData.obs.index,:]
    return RawData

def RunClustering(adata):
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca_variance_ratio(adata, log=True)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
#    sc.tl.tsne(adata)
    sc.tl.louvain(adata)
    sc.pl.umap(adata, color=['louvain',], legend_loc='on data')

def PlotPhen(adata, outpath, tp='phen'):
    sc.pl.umap(adata, color=['louvain'], alpha=0.6, size=12, save='ClassficationSep.pdf')
    sc.pl.umap(adata, color=['louvain'], alpha=0.6, size=12, legend_loc='on data',save='ClassficationOn.pdf' )
    sc.pl.umap(adata, color=[tp], alpha=0.6, size=12, save='ClassficationPhen.pdf')
    os.system('mv figures/umapClassficationSep.pdf {}'.format(outpath))
    os.system('mv figures/umapClassficationOn.pdf {}'.format(outpath))
    os.system('mv figures/umapClassficationPhen.pdf {}'.format(outpath))

def SelectPositiveCluster(adata, cutoff):
    obs = adata.obs
    dict_tmp = {}
    list_out = []
    for index in obs.index:
        louvain = obs.loc[index, 'louvain']
        if louvain not in dict_tmp:
            dict_tmp[louvain] = [obs.loc[index, 'phen']]
        else:
            dict_tmp[louvain] += [obs.loc[index, 'phen']]
    for key in sorted(dict_tmp.keys()):
        print(dict_tmp[key].count(1)/len(dict_tmp[key]))
        if dict_tmp[key].count(1)/len(dict_tmp[key]) >= cutoff:
            list_out.append(key)
    return list_out

def GetRankGenes(adata, list_cls, outpath, genes=20):
    sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon')
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    pd_out = pd.DataFrame(result['names']).loc[:str(genes-1),  list_cls]
    list_out = []
    for index in pd_out.index:
        for col in pd_out.columns:
            if pd_out.loc[index, col] not in list_out:
                if not pd.isnull(pd_out.loc[index, col]):
                    list_out.append(pd_out.loc[index, col])
    genes = open(os.path.join(outpath, 'FeatureGenes.lst'), 'w')
    for gene in sorted(list_out):
        genes.write(gene+'\n')
    genes.close()
    pd_x = adata[:, sorted(list_out)]
    pd_x = pd.DataFrame(pd_x.X, index=pd_x.obs.index, columns=pd_x.var.index, dtype='int')
    pd_x.to_csv(os.path.join(outpath, 'FeatureDiff.txt'), sep='\t', header=True, index=True)
    adata.write(os.path.join(outpath,'Feature.h5ad'))
    return pd_x

def main():
    parser = argparse.ArgumentParser(description="NG project Cluster pipeline")
    parser.add_argument('-m', help='input feature matrix, colomns is Gene ids, index is sample', required=True)
    parser.add_argument('-p', help='input phen matrix, colomns is sample phentype, index is sample', required=True)
    parser.add_argument('-c', help='cutoff', type=float, default=0.7)
    parser.add_argument('-o', help='output path', required=True)
    argv=vars(parser.parse_args())
    RawData = ReadData(argv['m'], argv['p'])
    RunClustering(RawData)
    PlotPhen(RawData, argv['o'])
    list_cls = SelectPositiveCluster(RawData, argv['c'])
    pd_x = GetRankGenes(RawData, list_cls, argv['o'])

if __name__ == '__main__':
    main()
