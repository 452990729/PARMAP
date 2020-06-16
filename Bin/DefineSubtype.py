#!/usr/bin/env python2


import sys
import re
import os
import argparse
import numpy as np
import pandas as pd
from HTSeq import FastaReader


def RefFa(file_in):
    dict_fa = {}
    list_all = []
    for item in FastaReader(file_in):
        list_lb = re.split('__', item.name)
        lb = list_lb[0]
        list_all.append(list_lb[1])
        dict_fa[lb] = str(item.seq)
    return sorted(set(list_all)), dict_fa

def MakeCluster(file_in, dict_fa, outpath):
    ClusterPath = os.path.join(outpath, 'Cluster')
    if os.path.exists(ClusterPath):
        os.system('rm -rf {}'.format(ClusterPath))
    os.mkdir(ClusterPath)
    m = 1
    with open(file_in, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                list_split = re.split('\t', line.strip('\n'))
                list_seq = [re.split('\(', list_split[0])[0],]+[re.split('\(', i)[0] for i in re.split('\s+', list_split[2])]
                outfile = open(os.path.join(ClusterPath, 'Cluster{}.fa'.format(str(m))), 'w')
                for seq in list_seq:
                    if len(seq) > 0:
                        outfile.write('>{}\n{}\n'.format(seq, dict_fa[seq]))
                outfile.close()
                m += 1

def MakeAllSample(list_in, list_all, label):
    list_out = []
    if '_' in list_in[0]:
        set_clip = set([re.split('_', m)[0] for m in list_in])
    elif 'fig' in list_in[0]:
        set_clip = set([re.findall('fig\|(\d+\.\d+)\.peg', m)[0] for m in list_in])
    for i in list_all:
        if i not in set_clip:
            list_out.append(0)
        else:
            list_out.append(1)
    return pd.DataFrame(np.array(list_out).reshape(1, len(list_out)),\
                        columns=list_all, index=[label+'_'+list_in[0],])

def HandleData(file_in, list_all, outfile):
    dict_fa = {}
    Cluster = os.path.basename(file_in).strip('.fa')
    for item in FastaReader(file_in):
        seq = str(item.seq)
        if seq not in dict_fa:
            dict_fa[seq] = [item.name,]
        else:
            dict_fa[seq] += [item.name,]
    n = 1
    for seq in dict_fa:
        pd_out = MakeAllSample(dict_fa[seq], list_all, Cluster)
        for index in pd_out.index:
            outfile.write(index+'\t'+'\t'.join([str(i) for i in pd_out.loc[index,:]])+'\n')

def SelectAll(path_in, list_all, outpath):
    outfile = open(os.path.join(outpath, 'tmp.txt'), 'w')
    outfile.write(' \t'+'\t'.join(list_all)+'\n')
    root, dirs, files = next(os.walk(path_in))
    for fl in files:
        HandleData(os.path.join(root, fl), list_all, outfile)
    outfile.close()
    pd_out = pd.read_csv(os.path.join(outpath, 'tmp.txt'), sep='\t', index_col=0, header=0)
    pd_out.T.to_csv(os.path.join(outpath, 'SubType.txt'), sep='\t', header=True, index=True)
    os.system('rm {}'.format(os.path.join(outpath, 'tmp.txt')))

def main():
    parser = argparse.ArgumentParser(description="Define the subtype of gene allele")
    parser.add_argument('-p', help='the input fasta pep file', required=True)
    parser.add_argument('-c', help='the input cluster file', required=True)
    parser.add_argument('-o', help='the output path', required=True)
    argv=vars(parser.parse_args())
    RawPep = argv['p']
    GeneCluster = argv['c']
    Outpath = argv['o']
    list_all, dict_fa = RefFa(RawPep)
    MakeCluster(GeneCluster, dict_fa, Outpath)
    SelectAll(os.path.join(Outpath, 'Cluster'), list_all, Outpath)


if __name__ == '__main__':
    main()




