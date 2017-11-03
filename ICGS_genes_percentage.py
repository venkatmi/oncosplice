#!/usr/bin/env python

import numpy as np
import pylab as pl
import sys,string
import os
import os.path
global test_list
import timeit
import time
from numpy import corrcoef, sum, log, arange
from scipy.stats.stats import pearsonr
import multiprocessing
import misopy as mi
from misopy import index_gff
import scipy.cluster.hierarchy as sch
from collections import defaultdict

if __name__ == '__main__':
    #manager=Manager()
    #start = timeit.default_timer()
    #correc=manager.dict()
   # create_sub_file()
   

    Splicecluster= defaultdict(list)
    countcluster=defaultdict(int)
    genelst=[]
    totalgene=[]
    splicelst=[]
    
    test_name="/Volumes/MyPassport/AIM1_MarkerGenes/Genes_groups.txt"
    text_file = open(test_name,'rU')
    lines = text_file.readlines()
    text_file.close()
    head=0
    for i in lines:
        if head==0:
            head=1
            continue
        i = i.rstrip(os.linesep)
       
        i=string.split(i,'\t')
        Splicecluster[i[1]].append(i[0])
        
    
    text_file1=open('/Volumes/MyPassport/Figure1_final/ExpressionInput/SamplePrediction/DataPlots/Clustering-Uncorrelatedevents_0.4-filtered1-CORRELATED-FEATURES-hierarchical_cosine_euclidean.txt','rU')
    lines1 = text_file1.readlines()
    head=0
    count=0
    for i in lines1:
        if head==0:
            head=1
            continue
        if 'column_clusters-flat' in i:
            continue
        i = i.rstrip(os.linesep)
        i=string.split(i,'\t')
       
        for key in Splicecluster:
            if i[0] in Splicecluster[key]:
               # print i[0]
                countcluster[key]=countcluster[key]+1
        count=count+1
    export_object2 = open("/Volumes/MyPassport/Figure1_final/ExpressionInput/SamplePrediction/DataPlots/Clustering-Uncorrelatedevents_0.4-filtered1-CORRELATED-FEATURES-hierarchical_cosine_euclidean-res.txt",'w')
    for key1 in countcluster:
        export_object2.write(key1+'\t'+str(countcluster[key1])+'\t'+str(len(Splicecluster[key1]))+'\n')
        
        
                
       
    
    
    
        
        