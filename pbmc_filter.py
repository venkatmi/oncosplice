#!/usr/bin/env python


import os
import os.path
global test_list
import matplotlib.pyplot as plt
import copy
from collections import Counter
from collections import defaultdict
import string
import numpy
cluster_grps = defaultdict(list)
import random
head=0
import math
samples_name=[]
samples=[]


def initiaterec():
    
    text_file = open("/Volumes/MyPassport/Singlecelldata/33k_CPTT_matrix-CORRELATED-FEATURES.txt", "rU")
    lines = text_file.readlines()
    text_file.close()
    head=0
    for line in lines:
        lin=line.rstrip('\r\n')
        t = string.split(lin,'\t')
        for i in t:
            samples.append(i)
    
        
    export=open("/Volumes/MyPassport/Singlecelldata/33k_CPTT_matrix-200_v3_upd.txt","w")
   # head=0
    count=1
    for line in open('/Volumes/MyPassport/Singlecelldata/33k_CPTT_matrix-200_v3.txt','rU').xreadlines():
        lin=line.rstrip('\r\n')
        t = string.split(lin,'\t')
        export.write(samples[count]+'\t'+str(t[0])+'\n')
        count=count+1
    #    if head==0:
    #        lin=line.rstrip('\r\n')
    #        t = string.split(lin,'\t')
    #        export.write
    #        for i in t:
    #            export.write('\t'+i)
    #        export.write('\n')
    #        head=1
    #        continue
    #    else:
    #        #if "HES4" in line:
    #         #   print line
    #            
    #       # lin=line.rstrip('\r\n')
    #        #t = string.split(lin,'\t')
    #       # samples.append(t[0])
    ##samples_names=random.sample(samples,15000)
    ##
    ##head=0
    ##        #print len(samples)
    ##for line in lines:
    ##
    ##    if head==0:
    ##        head=1
    ##        continue
    ##    lin=line.rstrip('\r\n')
    ##    t = string.split(lin,'\t')
    ##    
    ##    if str(t[0]) in samples_names:
    ##        #print t[0]
    #        export.write(line)
            
        #else:
        #    continue
           

initiaterec()
