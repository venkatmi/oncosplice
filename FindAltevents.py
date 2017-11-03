#!/usr/bin/env python


import numpy as np
import sys,string
import os
import os.path
from numpy import corrcoef, sum, log, arange
from scipy.stats.stats import pearsonr
import traceback
tandem = dict()
dem=dict()
lst=[]
count=0


genefile="/Users/meenakshi/Documents/ExpressionInput/PSIAnnotFile_upd_complete1.txt"
#altfile="/Users/meenakshi/Documents/ExpressionInput/PSIAnnotFile_upd_complete1_noalt.txt"




export=open("/Users/meenakshi/Documents/ExpressionInput/PSIAnnotFile_upd_complete1_sub.txt",'a')
head="true"
for line1 in open(genefile,'rU').xreadlines():
    newlist=[]
    line1 = line1.rstrip(os.linesep)
    if head=="true":
        head="false"
        continue
    t=string.split(line1,'\t')
    splice1=t[0]
    gene=t[1]
    iso1=t[6]
    iso2=t[7]
    count+=1
    if splice1 not in lst:
        export.write(splice1+'\n')
    flag="true"
    for line2 in open(genefile,'rU').xreadlines():
        if flag=="true":
            if splice1 in line2:
                flag="false"
            continue
        else:
            if t[1] in line2:
                t1=string.split(line2,'\t')
                splic11=t1[0]
                iso11=t1[6]
                iso12=t1[7]
                if iso1==iso11 or iso1==iso12:
                    lst.append(splic11)
                if iso2==iso11 or iso2==iso12:
                    lst.append(splic11)
            else:
                break
    if count%500==0:
        print count
    
    