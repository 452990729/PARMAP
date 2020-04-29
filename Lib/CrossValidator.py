#!/usr/bin/env python2


import re
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from scipy import interp, stats
from sklearn.preprocessing import label_binarize
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, auc, accuracy_score, f1_score, recall_score, precision_score
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
from compare_auc import roc_variance


def PlotROC(ax, lb, color, fpr, tpr, auc, ci):
    if len(ci)>0:
        ax.plot(fpr, tpr, 'k-', lw=2, label='{}AUC(%95 CI) = {}({})'.\
               format(lb, str(round(auc, 2)), '-'.join([str(round(i, 2)) for i in ci])), color=color)
    else:
        ax.plot(fpr, tpr, 'k-', lw=2, label='{}AUC = {}'.format(lb, str(round(auc, 2))), color=color)
    return ax

def PlotPrecisionRecall(lb, precision, recall, average_precision):
    plt.style.use(['my-paper'])
    fig, ax = plt.subplots(figsize=(10,8))
    ax.step(recall, precision, color='b', alpha=0.2, where='post')
    ax.fill_between(recall, precision, alpha=0.2, color='b')
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_ylim([0.0, 1.05])
    ax.set_xlim([0.0, 1.0])
    ax.set_title('Precision-Recall curve: AP={0:0.2f}'.format(
                  average_precision))
    fig.savefig('{}_PrecisionRecall.pdf'.format(lb))

def MultiClass(feature, response, model, fold):
    skf = StratifiedKFold(n_splits=fold, random_state=None, shuffle=False)
    response_befor = response
    response = pd.DataFrame(label_binarize(response, classes=list(set(list(response)))))
    n_classes = response.shape[1]
    micro_mean_tpr = 0.0
    micro_mean_fpr = np.linspace(0, 1, 100)
    macro_mean_tpr = 0.0
    macro_mean_fpr = np.linspace(0, 1, 100)
    for train_index, test_index in skf.split(feature, response_befor):
        fpr = {}
        tpr = {}
        print np.array(response.iloc[train_index])
        mbs = model.fit(feature.iloc[train_index], np.array(response.iloc[train_index]),)
        y_score = mbs.predict_proba(feature.iloc[test_index])[:,1]
        y_test = np.array(response.iloc[test_index])
        for i in range(n_classes):
            fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i])
        fpr["micro"], tpr["micro"], _  = roc_curve(y_test.ravel(), y_score.ravel())
        all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))
        mean_tpr = np.zeros_like(all_fpr)
        for i in range(n_classes):
            mean_tpr += interp(all_fpr, fpr[i], tpr[i])
        mean_tpr /= n_classes
        fpr["macro"] = all_fpr
        tpr["macro"] = mean_tpr
        micro_mean_tpr += interp(micro_mean_fpr, fpr["micro"], tpr["micro"])
        macro_mean_tpr += interp(micro_mean_fpr, fpr["macro"], tpr["macro"])
        micro_mean_tpr[0] = 0.0
        macro_mean_tpr[0] = 0.0
    micro_mean_tpr /= fold
    macro_mean_tpr /= fold
    micro_mean_tpr[-1] = 1.0
    macro_mean_tpr[-1] = 1.0
    micro_auc = auc(micro_mean_fpr, micro_mean_tpr)
    macro_auc = auc(macro_mean_fpr, macro_mean_tpr)
    return micro_mean_fpr, micro_mean_tpr, micro_auc,\
            macro_mean_fpr, macro_mean_tpr, macro_auc

def BinaClass(feature, response, model, fold):
    skf = StratifiedKFold(n_splits=fold, random_state=None, shuffle=False)
    list_y_score = []
    list_y_test = []
    list_y_pred = []
    list_accuracy = []
    list_precision = []
    list_recall = []
    list_f1 = []
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
#    response = response.dropna()
#    feature = feature.dropna()
#    feature = feature.loc[response.index, :]
#    response = response.loc[feature.index]
    for train_index, test_index in skf.split(feature, response):
        mbs = model.fit(feature.iloc[train_index], response.iloc[train_index],)
        y_score = mbs.predict_proba(feature.iloc[test_index])[:,1]
        y_test = response.iloc[test_index]
        y_pred = mbs.predict(feature.iloc[test_index])
        list_accuracy.append(accuracy_score(y_test, y_pred))
        list_precision.append(precision_score(y_test, y_pred))
        list_recall.append(recall_score(y_test, y_pred))
        list_f1.append(f1_score(y_test, y_pred))
        list_y_score += list(y_score)
        list_y_test += list(np.array(y_test.T))
        list_y_pred += list(y_pred)
        fpr, tpr, _ = roc_curve(y_test, y_score)
        mean_tpr += interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
    mean_tpr /= fold
    mean_tpr[-1] = 1.0
    mean_accuracy = round(sum(list_accuracy)/fold, 4)
    mean_precision = round(sum(list_precision)/fold, 4)
    mean_recall = round(sum(list_recall)/fold, 4)
    mean_f1 = round(sum(list_f1)/fold, 4)
    auc_value, ci = GetAUC(np.array(list_y_test), np.array(list_y_score))
    out = open('{}_Score.txt'.format(re.split('\(', str(model))[0]), 'w')
    out.write(' \tauc\taccuracy\tprecision\trecall\tf1\n')
    out.write('\t'.join([re.split('\(', str(model))[0], str(auc_value), str(mean_accuracy),
                         str(mean_precision), str(mean_recall), str(mean_f1)]))
    out.close()
    precision, recall, _ = precision_recall_curve(list_y_test, list_y_pred)
    average_precision = average_precision_score(list_y_test, list_y_pred)
    return mean_fpr, mean_tpr, auc_value, ci, precision, recall, average_precision

def GetAUC(np_test, np_predict):
    alpha=0.95
    auc_value, auc_cov = delong_roc_variance(np_test, np_predict)
    auc_std = np.sqrt(auc_cov)
    lower_upper_q = np.abs(np.array([0, 1]) - (1 - alpha) / 2)
    ci = stats.norm.ppf(lower_upper_q, loc=auc_value, scale=auc_std)
    return auc_value,ci

def Predict(feature, response, model):
    y_score = model.predict_proba(feature)[:,1]
    y_test = response
    y_pred = model.predict(feature)
    fpr, tpr, _ = roc_curve(y_test, y_score)
    accuracy_value = accuracy_score(y_test, y_pred)
    precision_value = precision_score(y_test, y_pred)
    recall_value = recall_score(y_test, y_pred)
    f1_value = f1_score(y_test, y_pred)
    precision, recall, _ = precision_recall_curve(y_test, y_pred)
    average_precision = average_precision_score(y_test, y_pred)
    auc_value, ci = GetAUC(np.array(list(np.array(y_test.T))), np.array(list(y_score)))
    out = open('{}_Score.txt'.format(re.split('\(', str(model))[0]), 'w')
    out.write(' \tauc\taccuracy\tprecision\trecall\tf1\n')
    out.write('\t'.join([re.split('\(', str(model))[0], str(auc_value), str(accuracy_value),
                         str(precision_value), str(recall_value), str(f1_value)]))
    out.close()
    return fpr, tpr, auc_value, ci, precision, recall, average_precision

def MakeROC(tp, feature, response, model, test_x=None, test_y=None, fold=5, bootstrap=1):
    plt.style.use(['my-paper', 'my-line'])
    colors = ['darkorange', 'blue', 'red', 'yellow']
    fig, ax = plt.subplots(figsize=(10,8))
    ax.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6))
    ax.set_xlim([-0.05, 1.05])
    ax.set_ylim([-0.05, 1.05])
    ax.set_xlabel('1-Specificity')
    ax.set_ylabel('Sensitivity')
    ax.spines['right'].set_visible(True)
    ax.spines['top'].set_visible(True)
    if tp in ['multi', 'bina']:
        m = 0
        if bootstrap == 1:
            for col in response.columns:
                RunValidator(tp, feature, response.loc[:,col], model, fold, ax, col, colors[m])
                m += 1
        else:
            RunValidator(tp, feature, response.iloc[:,0], model, fold, ax, '', colors[m], bootstrap)
            m += 1
    elif tp == 'predict':
        m = 0
        for col in response.columns:
            mbs = model.fit(feature, response.iloc[:,m],)
            fpr, tpr, auc_value, ci, precision, recall, average_precision = Predict(test_x, test_y.loc[:, col], mbs)
            PlotROC(ax, col, colors[m], fpr, tpr, auc_value, ci)
            PlotPrecisionRecall(re.split('\(', str(model))[0], precision, recall, average_precision)
            m += 1
    ax.legend(loc="lower right")
    fig.savefig('{}_ROC.pdf'.format(re.split('\(', str(model))[0]))

def RunValidator(tp, feature, response, model, fold, ax, lb, color, bootstrap=1):
    if tp == 'multi':
        micro_mean_fpr, micro_mean_tpr, micro_auc, macro_mean_fpr, macro_mean_tpr, macro_auc =\
                MultiClass(feature, response, model, fold=fold)
        PlotROC(ax, 'micro-average ', 'deeppink', micro_mean_fpr, micro_mean_tpr, micro_auc)
        PlotROC(ax, 'macro-average ', 'navy', macro_mean_fpr, macro_mean_tpr, macro_auc)
    elif tp == 'bina':
        if bootstrap == 1:
            mean_fpr, mean_tpr, mean_auc, ci, precision, recall, average_precision = BinaClass(feature, response, model, fold)
            PlotROC(ax, lb, color, mean_fpr, mean_tpr, mean_auc, ci)
            PlotPrecisionRecall(re.split('\(', str(model))[0], precision, recall, average_precision)
        else:
            for i in range(bootstrap):
                if i == 0:
                    mean_fpr, mean_tpr, mean_auc, ci = BinaClass(feature, response, model, fold)
                    auc = mean_auc
                    ci = ''
                    PlotROC(ax, lb, color, mean_fpr, mean_tpr, auc, ci)
                else:
                    mean_fpr, mean_tpr, mean_auc, ci = BinaClass(feature, response, model, fold)
                    auc = (auc+mean_auc)/2
                    ci = ''
                    PlotROC(ax, lb, color, mean_fpr, mean_tpr, auc, ci)
