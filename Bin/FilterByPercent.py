#!/usr/bin/env python2


import sys
import re
import argparse
import numpy as np
import pandas as pd

def ReadData(file_in):
    pd_data = pd.read_csv(file_in, sep='\t', header=0, index_col=0)
    return pd_data

def CallMix(value, pd_sus):
    total = 0
    for vls in pd_sus.values:
        total += (float(vls)/value)**2
    return round(total/pd_sus.shape[0], 2)

def MakeCut(pd_data, pd_cls, ctl):
    list_r = [str(i) for i in pd_cls.loc[pd_cls['Type']=='R', :].index]
    list_s = [str(i) for i in pd_cls.loc[pd_cls['Type']=='S', :].index]
    pd_r = pd_data.loc[:, list_r]
    pd_s = pd_data.loc[:, list_s]
    list_tmp = []
    for index in pd_r.index:
        value = pd_r.loc[index,:].max()
        pd_sus = pd_s.loc[index,:]
        if value >= 0.6:
            if CallMix(value, pd_sus) <= ctl:
                list_tmp.append(index)
    pd_out = pd_data.loc[list_tmp, :]
    return pd_out

def main():
    parser = argparse.ArgumentParser(description="filter feature by Percent")
    parser.add_argument('-m', help='input matrix data', required=True)
    parser.add_argument('-c', help='input class file', required=True)
    parser.add_argument('-t', help='cutoff<<0.1>>', type=float, default=0.1)
    parser.add_argument('-o', help='output file<<FeatureFilterByPercent.txt>>', default='FeatureFilterByPercent.txt')
    argv=vars(parser.parse_args())
    pd_data = ReadData(argv['m'])
    pd_cls = ReadData(argv['c'])
    pd_out = MakeCut(pd_data, pd_cls, argv['t'])
    pd_out.to_csv(argv['o'], sep='\t', header=True, index=True)


if __name__ == '__main__':
    main()
