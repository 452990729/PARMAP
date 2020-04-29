# PARMAP
A pan-genome based computational framework to predict antimicrobial resistance

| | |
|---|---|
| Author | Xuefei Lee |
| License | See included LICENSE |

Bacterial resistance is extremely serious, and multi-drug resistant bacteria are emerging in endlessly. The development of efficient computational biology tools is very important for the prediction of bacterial resistance. PARMAP is a computing framework for predicting bacterial resistance. It predicts the resistance of bacteria to multiple antibiotics based on the coding sequence of the strain's gene protein.

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
Step2:align the pangeome data.
```bash
./AlignCluster.py -h
usage: AlignCluster.py [-h] -c C -f F [-t T] -o O [-r]

Align The Corpan Cluster By Muscle

optional arguments:
  -h, --help  show this help message and exit
  -c C        the input cluster file
  -f F        the input fasta file
  -t T        the core used
  -o O        the output path
  -r          run now
```
Step3:get the gene allele matrix.
```bash
./DefineSubtype.py -h
usage: DefineSubtype.py [-h] -m M -c C -o O

define pan-genome Subtype

optional arguments:
  -h, --help  show this help message and exit
  -m M        input all sample pep fasta file ,obtained from pan genome
              analysis pipeline
  -c C        pan-genome cluster, obtained from pan genome analysis pipeline
  -o O        output path
```
Step4:filter gene allele matrix and get the key feature.
```bash
./FilterSubtype.py -h
usage: FilterSubtype.py [-h] -m M -p P [-c C] -o O

NG project Cluster pipeline

optional arguments:
  -h, --help  show this help message and exit
  -m M        input feature matrix, colomns is Gene ids, index is sample
  -p P        input phen matrix, colomns is sample phentype, index is sample
  -c C        cutoff
  -o O        output path
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


Citation
--------


License
-------
Comply with GNU GENERAL PUBLIC LICENSE.
