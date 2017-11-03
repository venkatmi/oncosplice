#!/usr/bin/env python

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
lis=[]

def strip_first_col(fname, delimiter=None):
    with open(fname, 'r') as fin:
        for line in fin:
            try:
               yield line.split(delimiter, 1)[1]
            except IndexError:
               continue

Xobs = np.loadtxt(strip_first_col("/Volumes/MyPassport/Leucegene/Leucegene_Results/AlternativeOutput/ExpressionInput/exp.LeucegeneSplicing-filtered.txt"),skiprows=1)
X=np.loadtxt("/Volumes/MyPassport/Leucegene/Leucegene_Results/AlternativeOutput/ExpressionInput/U2AF1like_SRSF2like_v2.txt")
#beta=zip(*X)
X=X.reshape(len(X), 1)
#beta=np.array(beta)
Y=np.array(Xobs)
Xobs=zip(*Xobs)
Xobs=np.array(Xobs)
export=open('test_v2.txt','w')


#lst=[]

#pc=decomposition.PCA(n_components=1)
#pc.fit(Y)
#a=pc.components_[[0]][0]*1000
#np.savetxt("pccoeff.txt",a)
#print len(a)
regr = linear_model.LinearRegression()
print X.shape,Xobs.shape
regr.fit(X,Xobs)
print('Coefficients: \n', regr.coef_)
np.savetxt("coefficients.txt",regr.coef_)
print regr.coef_.shape
regr.coef_=np.matrix(regr.coef_)
#X=np.matrix(X)
X=X.reshape(1,len(X))
k=regr.coef_*X
xlen,ylen= k.shape
Xobs=0;
nestedlst=[]
for i in range(xlen):
    lst=[]
    for j in range(ylen):
        if Y[i,j]==0:
            lst.append('')
        else:
            final=Y[i,j]-k[i,j]
            
            final=round(final,5)
            lst.append(str(final))
    nestedlst.append(lst)
    
        


export.writelines('\t'.join(i) + '\n' for i in nestedlst)

    
    
    





