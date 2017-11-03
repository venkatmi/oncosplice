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
from sklearn.decomposition import SparsePCA
lis=[]

def strip_first_col(fname, delimiter=None):
    with open(fname, 'r') as fin:
        for line in fin:
            try:
               yield line.split(delimiter, 1)[1]
            except IndexError:
               continue

X = np.loadtxt(strip_first_col("/Volumes/MyPassport/Leucegene/Leucegene_Results/AlternativeOutput/ExpressionInput/exp.LeucegeneSplicing.txt"),skiprows=1)
X=zip(*X)
X=np.array(X)
sparsepc=decomposition.SparsePCA()
sparsepc.fit(X)
print("Optimal number of features : %d" % sparsepc.components_)


