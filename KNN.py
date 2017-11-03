#!/usr/bin/env python
from sklearn.neighbors import KNeighborsClassifier
import numpy as np
import pylab as pl
import heapq
import sys,string
from numpy import loadtxt
import os
import os.path
global test_list
import matplotlib.pyplot as plt
from sklearn import decomposition
from sklearn.decomposition import PCA
from sklearn import datasets, linear_model
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import make_moons, make_circles, make_classification
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.ensemble import BaggingClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.lda import LDA
from sklearn.qda import QDA
from numpy import loadtxt
#from sklearn import cross_validation
from sklearn.cross_validation import StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn.metrics import zero_one_loss
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_score
from sklearn import svm
from sklearn.multiclass import OneVsOneClassifier
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import LinearSVC
from sklearn import linear_model
lis=[]

def strip_first_col(fname, delimiter=None):
    with open(fname, 'r') as fin:
        for line in fin:
            try:
               yield line.split(delimiter, 1)[1]
            except IndexError:
               continue
#Load Training
Xobs = np.loadtxt(strip_first_col("/Volumes/MyPassport/KNN_classification/ProperclassificationTest/GeneExpression/Classification/Train_Data.txt"),skiprows=1)
Xobs=zip(*Xobs)
#Load groups
X=np.loadtxt("/Volumes/MyPassport/KNN_classification/ProperclassificationTest/GeneExpression/Classification/groups.train.txt")
print len(Xobs)
#Load test
Y=np.loadtxt("/Volumes/MyPassport/KNN_classification/ProperclassificationTest/GeneExpression/Classification/Test_Data.txt")
Y=zip(*Y)
#Reshape the files
X=X.reshape(len(X), 1)
Xobs=np.array(Xobs)
Y=np.array(Y)


regr = KNeighborsClassifier(n_neighbors=1)
print X.shape,Xobs.shape
regr.fit(Xobs,X[:,0])
q=regr.predict(Y)
np.savetxt("/Volumes/MyPassport/KNN_classification/ProperclassificationTest/GeneExpression/Classification/complete_KNN.txt",q)

regr = LinearSVC()
regr.fit(Xobs,X[:,0])
q=regr.predict(Y)
np.savetxt("/Volumes/MyPassport/KNN_classification/ProperclassificationTest/GeneExpression/Classification/complete_SVC.txt",q)


regr = RandomForestClassifier(n_estimators=60,max_features='sqrt',bootstrap='true')
regr.fit(Xobs,X[:,0])
q=regr.predict(Y)
np.savetxt("/Volumes/MyPassport/KNN_classification/ProperclassificationTest/GeneExpression/Classification/complete_random_Forest.txt",q)

regr = LDA()
regr.fit(Xobs,X[:,0])
q=regr.predict(Y)
np.savetxt("/Volumes/MyPassport/KNN_classification/ProperclassificationTest/GeneExpression/Classification/complete_LDA.txt",q)      

#regr = LinearSVC()
#regr.fit(Xobs,X[:,0])
#q=regr._predict_proba_lr(Y)
#qr=regr.predict(Y)
#np.savetxt("/Users/meenakshi/Documents/KNN_classification/GeneExpression/results/zrsr2_proba_lr_linearSVC.txt",q)
#np.savetxt("/Users/meenakshi/Documents/KNN_classification/results/complete_predict_linearSVC.txt",qr)
    
        




