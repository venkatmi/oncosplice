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



genefile="/Volumes/MyPassport/U2AF1likeSplicingFactorTranscriptioncorrelations/splicingdict1.txt"
corrfile="/Volumes/MyPassport/Leucegene_U2AF1vsU2AF1like/result_binarized_0.8_pos.txt"

for line in open(genefile,'rU').xreadlines():
    line = line.rstrip(os.linesep)
  #  line=line.replace('\"','')
    t=string.split(line,'\t')
    
    gene=t[0]
    val=t[1]
    #print gene,val;sys.exit()
    tandem[gene]=val


export=open("correlation_0.8_pos_upd.txt",'w')

for line1 in open(corrfile,'rU').xreadlines():
    newlist=[]
    line1 = line1.rstrip(os.linesep)
    line1=line1.replace('\"','')
    t=string.split(line1,'\t')
    splice1=t[1]
    t[2]=t[2].replace("-",".")
    t[2]=t[2].replace("&",".")
    splice2=t[2]
    try:
    #print splice1,splice2;sys.exit()
        val1=tandem[splice1]
        val2=tandem[splice2]
        export.write(val1+'\t'+val2+'\t'+t[3]+'\n')
    except Exception:
        continue
    
    