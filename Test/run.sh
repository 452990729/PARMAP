#!/bin/bash

SCRIPT_DIR=$(cd $(dirname ${BASH_SOURCE[0]}); pwd)
pep=$SCRIPT_DIR/RawData/pep.lst
outpath=$SCRIPT_DIR/run

#===============================================================
mkdir $outpath
ln -s $SCRIPT_DIR/RawData/Phen 0.phen
$SCRIPT_DIR/../Bin/CorPanPipeline.py -c $pep -o $outpath
#===============================================================

mkdir 3.Align
cd 3.Align
$SCRIPT_DIR/../Bin/AlignCluster.py -c ../2.Stat/Sample_GeneCluster.list -f ../1.Cluster/all_gene.fa -o ./ -t 8 -r
$SCRIPT_DIR/../Bin/DefineSubtype.py ./Cluster ../0.phen/All SubType.txt
cd ..
#===============================================================

mkdir 4.Train
cd 4.Train
$SCRIPT_DIR/../Bin/ExtractById.py -m ../3.Align/SubType.txt -i ../0.phen/Train -o Feature.txt

$SCRIPT_DIR/../Bin/FilterSubtype.py -m Feature.txt -p ../0.phen/Train -o ./
mkdir ML
cd ML
$SCRIPT_DIR/../Bin/MlAndPredict.py validator -x ../FeatureDiff.txt -y ../../0.phen/Train
cd ../../
#===============================================================

mkdir 5.Test
cd 5.Test
$SCRIPT_DIR/../Bin/ExtractById.py -m ../3.Align/SubType.txt -i ../0.phen/Test -o Feature.txt
$SCRIPT_DIR/../Bin/ExtractById.py -m Feature.txt -i ../4.Train/FeatureGenes.lst -o FeatureDiff.txt -col
mkdir ML
cd ML
$SCRIPT_DIR/../Bin/MlAndPredict.py predict -x ../../4.Train/FeatureDiff.txt -y ../../0.phen/Train -testx ../FeatureDiff.txt -testy ../../0.phen/Test -m logistic
$SCRIPT_DIR/../Bin/MlAndPredict.py predict -x ../../4.Train/FeatureDiff.txt -y ../../0.phen/Train -testx ../FeatureDiff.txt -testy ../../0.phen/Test -m svm
$SCRIPT_DIR/../Bin/MlAndPredict.py predict -x ../../4.Train/FeatureDiff.txt -y ../../0.phen/Train -testx ../FeatureDiff.txt -testy ../../0.phen/Test -m rf
$SCRIPT_DIR/../Bin/MlAndPredict.py predict -x ../../4.Train/FeatureDiff.txt -y ../../0.phen/Train -testx ../FeatureDiff.txt -testy ../../0.phen/Test -m gdbt
cd ../../
#===============================================================

