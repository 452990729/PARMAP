#!/usr/bin/env python3
# coding: utf-8


import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import re
import os
import sys
import argparse
from matplotlib import pyplot as plt
import seaborn as sns
import scipy as sp
from scipy import stats
from scipy import interpolate
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=300, dpi_save=600, vector_friendly=False)

def ReadData(feature, phen):
    RawData = ad.read_csv(feature, delimiter='\t', first_column_names=True)
    phen = pd.read_csv(phen, sep='\t', header=0, index_col=0)
    phen.index = [str(i) for i in phen.index]
    RawData.obs['phen'] = phen.loc[RawData.obs.index,:]
    return RawData

def RunClustering(adata):
    sc.tl.pca(adata, svd_solver='arpack')
#    sc.pl.pca_variance_ratio(adata)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
#    sc.tl.tsne(adata)
    sc.tl.louvain(adata)
#    sc.pl.umap(adata, color=['louvain',], legend_loc='on data')

def PlotPhen(adata, outpath='figures', tp='phen'):
    matplotlib.use('Agg')
    sc.pl.umap(adata, color=['louvain'], alpha=0.6, size=12, save='ClassficationSep.pdf')
    sc.pl.umap(adata, color=['louvain'], alpha=0.6, size=12, legend_fontsize=8, legend_loc='on data',save='ClassficationOn.pdf' )
    sc.pl.umap(adata, color=[tp], alpha=0.6, size=12, save='ClassficationPhen.svg')
    os.system('mv figures/umapClassficationSep.pdf {}'.format(outpath))
    os.system('mv figures/umapClassficationOn.pdf {}'.format(outpath))
    os.system('mv figures/umapClassficationPhen.svg {}'.format(outpath))

def SelectPositiveCluster(adata, cutoff=0.7):
    obs = adata.obs
    dict_tmp = {}
    list_out = []
    for index in obs.index:
        louvain = obs.loc[index, 'louvain']
        if louvain not in dict_tmp:
            dict_tmp[louvain] = [obs.loc[index, 'phen']]
        else:
            dict_tmp[louvain] += [obs.loc[index, 'phen']]
    for key in sorted([int(i) for i in dict_tmp.keys()]):
        key = str(key)
        pert = dict_tmp[key].count(1)/len(dict_tmp[key])
        print(key+'\t'+str(pert)+'\t'+str(1-pert))
        if dict_tmp[key].count(1)/len(dict_tmp[key]) >= cutoff:
            list_out.append(key)
    return list_out

def GetClusterPercent(adata, list_gene, ctl):
    if ctl:
        list_cls = list(set(list(adata.obs['louvain'])))
    else:
        list_cls = SelectPositiveCluster(adata)
    pd_out = pd.DataFrame(index=list_gene, columns=list_cls)
    for index in pd_out.index:
        for columns in pd_out.columns:
            list_sm = list(adata[adata.obs['louvain']==columns][:, index].X)
            pd_out.loc[index, columns] = round(float(list_sm.count(1))/len(list_sm), 2)
    return pd_out

def GetCluster(adata, list_p, outfile):
    list_tmp = []
    with open(outfile, 'w') as f:
        f.write(' \tType\n')
        for cls in list_p:
            f.write(cls+'\tR\n')
            list_tmp.append(cls)
        for cls in sorted(set(list(adata.obs['louvain']))):
            if cls not in list_p:
                f.write(cls+'\tS\n')
                list_tmp.append(cls)
    return list_tmp

def FisherTest(adata_cls, adata_noncls):
    pd_cls = pd.DataFrame(adata_cls.X, index=adata_cls.obs.index, columns=adata_cls.var.index)
    pd_noncls = pd.DataFrame(adata_noncls.X, index=adata_noncls.obs.index, columns=adata_noncls.var.index)
    pd_r = pd.DataFrame(index=pd_cls.columns)
    list_p = []
    for i in range(pd_cls.shape[1]):
        list_s1 = list(np.array(pd_cls.iloc[:, i].T))
        list_s2 = list(np.array(pd_noncls.iloc[:, i].T))
        a = list_s1.count(0)
        b = list_s2.count(0)
        c = list_s1.count(1)
        d = list_s2.count(1)
        table = np.array([[a,b],[c,d]])
        s, p = stats.fisher_exact(table, alternative='less')
        list_p.append(p)
    pd_r['pvalue']=list_p
    return pd_r

def estimate(pv, m=None, verbose=False, lowmem=False, pi0=None):
    """
    Estimates q-values from p-values

    Args
    =====

    m: number of tests. If not specified m = pv.size
    verbose: print verbose messages? (default False)
    lowmem: use memory-efficient in-place algorithm
    pi0: if None, it's estimated as suggested in Storey and Tibshirani, 2003.
         For most GWAS this is not necessary, since pi0 is extremely likely to be
         1

    """
    assert(pv.min() >= 0 and pv.max() <= 1), "p-values should be between 0 and 1"

    original_shape = pv.shape
    pv = pv.ravel()  # flattens the array in place, more efficient than flatten()

    if m is None:
        m = float(len(pv))
    else:
        # the user has supplied an m
        m *= 1.0

    # if the number of hypotheses is small, just set pi0 to 1
    if len(pv) < 100 and pi0 is None:
        pi0 = 1.0
    elif pi0 is not None:
        pi0 = pi0
    else:
        # evaluate pi0 for different lambdas
        pi0 = []
        lam = sp.arange(0, 0.90, 0.01)
        counts = sp.array([(pv > i).sum() for i in sp.arange(0, 0.9, 0.01)])
        for l in range(len(lam)):
            pi0.append(counts[l]/(m*(1-lam[l])))

        pi0 = sp.array(pi0)

        # fit natural cubic spline
        tck = interpolate.splrep(lam, pi0, k=3)
        pi0 = interpolate.splev(lam[-1], tck)
        if verbose:
            print("qvalues pi0=%.3f, estimated proportion of null features " % pi0)

        if pi0 > 1:
            if verbose:
                print("got pi0 > 1 (%.3f) while estimating qvalues, setting it to 1" % pi0)
            pi0 = 1.0

    assert(pi0 >= 0 and pi0 <= 1), "pi0 is not between 0 and 1: %f" % pi0

    if lowmem:
        # low memory version, only uses 1 pv and 1 qv matrices
        qv = sp.zeros((len(pv),))
        last_pv = pv.argmax()
        qv[last_pv] = (pi0*pv[last_pv]*m)/float(m)
        pv[last_pv] = -sp.inf
        prev_qv = last_pv
        for i in range(int(len(pv))-2, -1, -1):
            cur_max = pv.argmax()
            qv_i = (pi0*m*pv[cur_max]/float(i+1))
            pv[cur_max] = -sp.inf
            qv_i1 = prev_qv
            qv[cur_max] = min(qv_i, qv_i1)
            prev_qv = qv[cur_max]

    else:
        p_ordered = sp.argsort(pv)
        pv = pv[p_ordered]
        qv = pi0 * m/len(pv) * pv
        qv[-1] = min(qv[-1], 1.0)

        for i in range(len(pv)-2, -1, -1):
            qv[i] = min(pi0*m*pv[i]/(i+1.0), qv[i+1])

        # reorder qvalues
        qv_temp = qv.copy()
        qv = sp.zeros_like(qv)
        qv[p_ordered] = qv_temp

    # reshape qvalues
    qv = qv.reshape(original_shape)

    return qv

def GetTestResult(adata, list_cls, outpath, qv=0.05):
    set_all = set([])
    for cls in list_cls:
        adata_sub = adata[adata.obs['louvain']==cls,:]
        adata_noncls = adata[adata.obs['louvain']!=cls,:]
        pd_r = FisherTest(adata_sub, adata_noncls)
        q_value = estimate(np.array(pd_r.iloc[:,0].T))
        pd_r.insert(0, 'q_value', q_value)
        pd_filter = pd_r.loc[pd_r['q_value']<qv,:]
#        pd_filter = pd_r.loc[pd_r['p_value']<qv,:]
        set_all = set_all|set(pd_filter.index)
    adata_x = adata[:, sorted(list(set_all))]
    pd_x = pd.DataFrame(adata_x.X, index=adata_x.obs.index, columns=adata_x.var.index, dtype='int')
    pd_x.to_csv(os.path.join(outpath, 'FeatureDiff.txt'), sep='\t', header=True, index=True)
    adata.write(os.path.join(outpath,'Feature.h5ad'))
    genes = open(os.path.join(outpath, 'FeatureGenes.lst'), 'w')
    for gene in sorted(list(set_all)):
        genes.write(gene+'\n')
    genes.close()
    return pd_x

def main():
    parser = argparse.ArgumentParser(description="NG project Cluster pipeline")
    parser.add_argument('-m', help='input feature matrix, colomns is Gene ids, index is sample', required=True)
    parser.add_argument('-p', help='input phen matrix, colomns is sample phentype, index is sample', required=True)
    parser.add_argument('-o', help='output path', required=True)
    argv=vars(parser.parse_args())
    RawData = ReadData(argv['m'], argv['p'])
    RunClustering(RawData)
    PlotPhen(RawData, argv['o'])
    list_cls = SelectPositiveCluster(RawData)
    pd_makers = GetTestResult(RawData, list_cls, argv['o'])
    with open(os.path.join(argv['o'], 'FeatureGenes.lst'), 'r') as f:
        list_gene = [i.strip() for i in f]
    pd_percent_all = GetClusterPercent(RawData, list_gene, 1)
    list_sort_cls = GetCluster(RawData, list_cls, os.path.join(argv['o'], 'Cluster.txt'))
    pd_percent_all = pd_percent_all.loc[:, list_sort_cls]
    pd_percent_all.to_csv(os.path.join(argv['o'],'ClusterDiffGenePercent_all.txt'), sep='\t', header=True, index=True)


if __name__ == '__main__':
	main()
