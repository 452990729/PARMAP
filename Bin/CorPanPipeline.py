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

def HandleConfig(path_in):
    if path_in.startswith('..'):
        return os.path.join(BasePath, path_in)
    else:
        return path_in

#### SOFT
PYTHON = config.get('SOFTWARE', 'python')
PERL = config.get('SOFTWARE', 'perl')
CDHIT = config.get('SOFTWARE', 'cd-hit')

### SCRIPT
CorPan_get_all = config.get('SCRIPT', 'CorPan_get_all')
CorPan_stat_list = config.get('SCRIPT', 'CorPan_stat_list')
CorPan_gene = config.get('SCRIPT', 'CorPan_gene')
CorPan_specific = config.get('SCRIPT', 'CorPan_specific')

#### DATABASE

def WriteShell(in_file, outpath):
    ClusterPath = os.path.join(outpath, '1.Cluster')
    out_cluster_shell = os.path.join(ClusterPath, 'run.sh')
    StatPath = os.path.join(outpath, '2.Stat')
    out_stat_shell = os.path.join(StatPath, 'run.sh')
    if os.path.exists(ClusterPath):
        os.system('rm -rf {}'.format(ClusterPath))
    os.mkdir(ClusterPath)
    out_cluster = open(out_cluster_shell, 'w')
    out_cluster.write('cd {}\n'.format(ClusterPath))
    out_cluster.write('{} {} {} all_gene.fa 0\n'.format(PERL, CorPan_get_all, in_file))
    out_cluster.write('{} -i all_gene.fa -o Sample -M 0 -c 0.5 -n 3 -p 1 -T 4 -g 1 -d 0 -s 0.7 -aL 0.7 -aS 0.7\n'.\
                     format(CDHIT))
    out_cluster.write('mv Sample Sample.fa')
    out_cluster.close()
    if os.path.exists(StatPath):
        os.system('rm -rf {}'.format(StatPath))
    os.mkdir(StatPath)
    out_stat = open(out_stat_shell, 'w')
    out_stat.write('cd {}\n'.format(StatPath))
    out_stat.write('{} {} ../1.Cluster/Sample.clstr {} --prefix Sample --outdir . --min_identity 0.5 --min_Rcov 0.7 --min_Scov 0.7\n'.\
                  format(PERL, CorPan_stat_list, in_file))
    out_stat.write('{} {} Sample_Core.matrix ../1.Cluster/all_gene.fa Sample_Core.fa\n'.format(PERL, CorPan_gene))
    out_stat.write('{} {} Sample_Pan.matrix ../1.Cluster/all_gene.fa Sample_Pan.fa\n'.format(PERL, CorPan_gene))
    out_stat.write('{} {} Sample_Dispensable.matrix ../1.Cluster/all_gene.fa Sample_Dispensable.fa\n'.format(PERL, CorPan_gene))
    out_stat.write('{} {} Sample_Specific.list {} Sample_Specific'.format(PERL, CorPan_specific, in_file))
    out_stat.close()
    os.system('sh {} 1>{}/run.log 2>{}/run.err'.format(out_cluster_shell, ClusterPath, ClusterPath))
    os.system('sh {} 1>{}/run.log 2>{}/run.err'.format(out_stat_shell, StatPath, StatPath))


def main():
    parser = argparse.ArgumentParser(description="Micro CorPan pipeline")
    parser.add_argument('-c', help='the input fasta list abs path sample\\tpath\\n', required=True)
    parser.add_argument('-o', help='the output abs path', required=True)
    argv=vars(parser.parse_args())
    WriteShell(argv['c'], argv['o'])


if __name__ == '__main__':
    main()
