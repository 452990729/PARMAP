# PARMAP
A pan-genome based computational framework to predict antimicrobial resistance

| | |
|---|---|
| Author | Xuefei Lee |
| License | See included LICENSE |

Bacterial resistance is extremely serious, and multi-drug resistant bacteria are emerging in endlessly. The development of efficient computational biology tools is very important for the prediction of bacterial resistance. PARMAP is a computing framework for predicting bacterial resistance. It predicts the resistance of bacteria to multiple antibiotics based on the coding sequence of the strain.

Installation
------------
Simply pull from github:

```bash
git clone https://github.com/452990729/PARMAP.git
```


Usage
-----
Step1:run pangenome anlysis based on the protein sequence.

```bash
./CorPanPipeline.py -h
usage: CorPanPipeline.py [-h] -c C -o O

Micro CorPan pipeline

optional arguments:
  -h, --help  show this help message and exit
  -c C        the input fasta list abs path sample\tpath\n
  -o O        the output abs path
```
Step2:get the gene allele matrix.
```bash
./DefineSubtype.py -h
usage: DefineSubtype.py [-h] -p P -c C -o O

define pan-genome Subtype

optional arguments:
  -h, --help  show this help message and exit
  -p P        input all sample pep fasta file ,obtained from pan genome
              analysis pipeline
  -c C        pan-genome cluster, obtained from pan genome analysis pipeline
  -o O        output path
```
Step3:cluster sample and define S/R clusters, also filter feature by fisher testingfilter gene allele matrix and get the key feature.
```bash
./FilterSubtype.py -h
usage: FilterSubtype.py [-h] -m M -p P -o O

NG project Cluster pipeline

optional arguments:
  -h, --help  show this help message and exit
  -m M        input feature matrix, colomns is Gene ids, index is sample
  -p P        input phen matrix, colomns is sample phentype, index is sample
  -o O        output path
```
Step4 filter feature by AMRS and get the key feature.
```bash
./FilterByPercent.py -h
usage: FilterByPercent.py [-h] -m M -c C [-t T] [-o O]

filter feature by Percent

optional arguments:
  -h, --help  show this help message and exit
  -m M        input matrix data
  -c C        input class file
  -t T        cutoff<<0.1>>
  -o O        output file prefix<<FinalFeature>>
```
Step5:mathine learning process, which you can choose the best model and parameters based on cross-validation and then test the model performance.
```bash
./MlAndPredict.py -h
usage: MlAndPredict.py [-h] -x X -y Y [-testx TESTX] [-testy TESTY]
                       [-m {svm,rf,logistic,gdbt}] [-t {bina,multi}] [-f F]
                       [-b B]
                       {predict,validator}

Mathine learning Validator or Predict

positional arguments:
  {predict,validator}   the method used predict/validator

optional arguments:
  -h, --help            show this help message and exit
  -x X                  input feature matrix, col are features, index are
                        samples
  -y Y                  input response matrix, col is response, index are
                        samples
  -testx TESTX          input test feature matrix, col are features, index are
                        samples
  -testy TESTY          input test response matrix, col are response(allow
                        multi response), index are samples
  -m {svm,rf,logistic,gdbt}
                        if use predict mode, the method used <<svm>>
  -t {bina,multi}       if use validator mode, the type of data <<bina>>
  -f F                  if use validator mode, cross time <<5>>
  -b B                  if use validator mode, bootstrap time <<1>>
```
You can get more details from the run.sh script from test directory.

Citation
--------


License
-------
Comply with GNU GENERAL PUBLIC LICENSE.
