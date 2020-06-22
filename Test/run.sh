#!/bin/bash

workdir=$(cd $(dirname $0); pwd)
SCRIPT_DIR=$workdir
pep=$workdir/RawData/pep.lst
outpath=$workdir/run
#===============================================================
mkdir $outpath
cd $outpath
ln -s $workdir/RawData/Phen 0.phen
$SCRIPT_DIR/../Bin/CorPanPipeline.py -c $pep -o $outpath
#===============================================================

mkdir 3.Align
cd 3.Align
$SCRIPT_DIR/../Bin/DefineSubtype.py -p ../1.Cluster/all_gene.fa -c ../2.Stat/Sample_GeneCluster.list -o ./
cd ..
#===============================================================
mkdir 4.Train
cd 4.Train
$SCRIPT_DIR/../Bin/ExtractById.py -m ../3.Align/SubType.txt -i ../0.phen/Train -o Feature.txt

$SCRIPT_DIR/../Bin/FilterSubtype.py -m Feature.txt -p ../0.phen/Train -o ./
$SCRIPT_DIR/../Bin/FilterByPercent.py -m ClusterDiffGenePercent_all.txt -c Cluster.txt -t 0.1
$SCRIPT_DIR/../Bin/ExtractById.py -m FeatureDiff.txt -i FinalFeature.lst -col -o FinalFeature.txt
mkdir ML
cd ML
$SCRIPT_DIR/../Bin/MlAndPredict.py validator -x ../FinalFeature.txt -y ../../0.phen/Train
cd ../../
#===============================================================

mkdir 5.Test
cd 5.Test
$SCRIPT_DIR/../Bin/ExtractById.py -m ../3.Align/SubType.txt -i ../0.phen/Test -o Feature.txt
$SCRIPT_DIR/../Bin/ExtractById.py -m Feature.txt -i ../4.Train/FinalFeature.lst -o FeatureDiff.txt -col
mkdir ML
cd ML
$SCRIPT_DIR/../Bin/MlAndPredict.py predict -x ../../4.Train/FinalFeature.txt -y ../../0.phen/Train -testx ../FeatureDiff.txt -testy ../../0.phen/Test -m logistic
$SCRIPT_DIR/../Bin/MlAndPredict.py predict -x ../../4.Train/FinalFeature.txt -y ../../0.phen/Train -testx ../FeatureDiff.txt -testy ../../0.phen/Test -m svm
$SCRIPT_DIR/../Bin/MlAndPredict.py predict -x ../../4.Train/FinalFeature.txt -y ../../0.phen/Train -testx ../FeatureDiff.txt -testy ../../0.phen/Test -m rf
$SCRIPT_DIR/../Bin/MlAndPredict.py predict -x ../../4.Train/FinalFeature.txt -y ../../0.phen/Train -testx ../FeatureDiff.txt -testy ../../0.phen/Test -m gdbt
cd ../../
#===============================================================

