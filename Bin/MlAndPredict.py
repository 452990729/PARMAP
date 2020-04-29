#!/usr/bin/env python2


import sys
import os
import argparse
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegressionCV, LassoCV
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier

BasePath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(BasePath+'/../Lib/')
from CrossValidator import MakeROC

def GetData(file_in):
    pd_data =pd.read_csv(file_in, sep='\t', header=0, index_col=0)
    return pd_data

def RunLasso(feature, response):
    clf = LassoCV(alphas = [1, 0.1, 0.001, 0.0005]).fit(feature, response)
    coef = pd.Series(clf.coef_, index = feature.columns)
    coef.to_csv('Lasso.coef.txt', sep='\t', index=True, header=False)

def GetSVM():
    SVM = svm.SVC(C=0.5, probability=True, class_weight='balanced')
    return SVM

def GetRF():
    RandomForest = RandomForestClassifier(n_estimators=50, max_depth=10, oob_score=True, class_weight='balanced')
    return RandomForest

def GetGDBT():
    GDBT = GradientBoostingClassifier(n_estimators=100, learning_rate=0.1, max_depth=3)
    return GDBT

def GetLogistic():
    Logistic = LogisticRegressionCV(multi_class="ovr",fit_intercept=True,\
                               Cs=np.logspace(-2,2,20),cv=2,penalty="l2",\
                               solver="lbfgs",tol=0.01,class_weight='balanced')
    return Logistic

def MakePredict(train_x, train_y, test_x, test_y, model):
    MakeROC('predict', test_x, test_y, model, test_x=test_x, test_y=test_y)

def Validator(tp, x, y, model, fold, bootstrap):
    MakeROC(tp, x, y, model, fold=fold, bootstrap=bootstrap)

def main():
    parser = argparse.ArgumentParser(description="Mathine learning Validator or Predict")
    parser.add_argument('method', help='the method used predict/validator', choices=['predict', 'validator'], nargs=1)
    parser.add_argument('-x', help='input feature matrix, col are features, index are samples', required=True)
    parser.add_argument('-y', help='input response matrix, col is response, index are samples', required=True)
    parser.add_argument('-testx', help='input test feature matrix, col are features, index are samples')
    parser.add_argument('-testy', help='input test response matrix, col are response(allow multi response), index are samples')
    parser.add_argument('-m', help='if use predict mode, the method used <<svm>>', choices=['svm', 'rf', 'logistic', 'gdbt'], default='svm')
    parser.add_argument('-t', help='if use validator mode, the type of data <<bina>>', choices=['bina', 'multi'], default='bina')
    parser.add_argument('-f', help='if use validator mode, cross time <<5>>', type=int, default=5)
    parser.add_argument('-b', help='if use validator mode, bootstrap time <<1>>', type=int, default=1)
    argv=vars(parser.parse_args())
    feature = GetData(argv['x'])
    response = GetData(argv['y'])
    response = response.loc[feature.index,:]
    if argv['method'][0] == 'predict':
        testx = GetData(argv['testx'])
        testy = GetData(argv['testy'])
        testy = testy.loc[testx.index,:]
        if argv['m'] == 'svm':
            model = GetSVM()
        elif argv['m'] == 'rf':
            model = GetRF()
        elif argv['m'] == 'logistic':
            model = GetLogistic()
        elif argv['m'] == 'gdbt':
            model = GetGDBT()
        MakePredict(feature, response, testx, testy, model)
    elif argv['method'][0] == 'validator':
        Validator(argv['t'], feature, response, GetSVM(), argv['f'], argv['b'])
        Validator(argv['t'], feature, response, GetRF(), argv['f'], argv['b'])
        Validator(argv['t'], feature, response, GetLogistic(), argv['f'], argv['b'])
        Validator(argv['t'], feature, response, GetGDBT(), argv['f'], argv['b'])


if __name__ == '__main__':
    main()
