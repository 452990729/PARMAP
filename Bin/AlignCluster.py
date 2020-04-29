#!/usr/bin/env python2

import sys
import re
import os
import argparse
import ConfigParser
from glob import glob

BasePath = os.path.split(os.path.realpath(__file__))[0]
config = ConfigParser.ConfigParser()
config.read(BasePath+'/config.ini')

#### SOFT
MUSCLE = config.get('SOFTWARE', 'muscle')

### SCRIPT
CutFile = config.get('SCRIPT', 'CutFile')


def ReadFasta(file_in):
    dict_fa = {}
    m = 0
    with open(file_in, 'r') as f:
        for line in f:
            line = line.strip()
            if m%2 == 0:
                label = re.split('__', line.lstrip('>'))[0]
            elif m%2 == 1:
                dict_fa[label] = line
            m += 1
    return dict_fa

def ExtractCluster(file_in, outpath, dict_fa):
    ClusterPath = os.path.join(outpath, 'Cluster')
    if os.path.exists(ClusterPath):
        os.system('rm -rf {}'.format(ClusterPath))
    os.mkdir(ClusterPath)
    m = 1
    Shell = open(os.path.join(outpath, 'Run.sh'), 'w')
    with open(file_in, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                list_split = re.split('\t', line.strip('\n'))
                list_seq = [re.split('\(', list_split[0])[0],]+[re.split('\(', i)[0] for i in re.split('\s+', list_split[2])]
                PanClusterPath = os.path.join(ClusterPath, 'Cluster{}'.format(str(m)))
                os.mkdir(PanClusterPath)
                fasta = open(os.path.join(PanClusterPath, 'input.fa'), 'w')
                for seq in list_seq:
                    if len(seq) > 0:
                        fasta.write('>{}\n{}\n'.format(seq, dict_fa[seq]))
                fasta.close()
                Shell.write('{} -in {} -out {}\n'.format(MUSCLE, os.path.join(PanClusterPath, 'input.fa'),\
                                                      os.path.join(PanClusterPath, 'output.fa')))
                m += 1
    Shell.close()

def main():
    parser = argparse.ArgumentParser(description="Align The Corpan Cluster By Muscle")
    parser.add_argument('-c', help='the input cluster file', required=True)
    parser.add_argument('-f', help='the input fasta file', required=True)
    parser.add_argument('-t', help='the core used', default=2)
    parser.add_argument('-o', help='the output path', required=True)
    parser.add_argument('-r', help='run now',  action='store_true')
    argv=vars(parser.parse_args())
    dict_fa = ReadFasta(argv['f'])
    ExtractCluster(argv['c'], argv['o'], dict_fa)
    if argv['r']:
        os.system('nohup {} --core {} --shell {} --outpath {}&'.\
                 format(CutFile, int(argv['t']), os.path.join(argv['o'], 'Run.sh'),\
                       argv['o']))


if __name__ == '__main__':
    main()
