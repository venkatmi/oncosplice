#!/usr/bin/env python

import numpy as np
import sys,string
import os
import os.path
from numpy import corrcoef, sum, log, arange
from scipy.stats.stats import pearsonr
import traceback
from collections import Counter

inFile = '/Volumes/MyPassport/AML_Images_Dec21/Figure3/Splicingevents_allmutations_uncor.txt'
export=open('/Volumes/MyPassport/AML_Images_Dec21/Figure3/Splicingevents_fold_sign_uncor.txt','w')
head=1;head1=1
countgrp=Counter()
tandem=Counter()
grplist=[]
query_data = open(inFile,'rU')
lines = query_data.readlines()

for line in lines:
    if head==1:
        head=0
    
        continue
    line = line.rstrip('\r\n')
    psi=string.split(line,'\t')
    countgrp[psi[1]]+=1
    if psi[1] not in grplist:
        grplist.append(psi[1])
    for line1 in lines:
        if head==1:
            head=0
            
            continue
        line1 = line1.rstrip('\r\n')
        psi1=string.split(line1,'\t')
        if psi1[0]==psi[0] and psi[1]!=psi1[1]:
            if (float(psi[2])>0.0 and float(psi1[2])>0.0) or (float(psi[2])<0.0 and float(psi1[2])<0.0):
                tandem[psi[1],psi1[1]]+=1

for i in grplist:
    export.write('\t'+i)
export.write('\n')

print tandem['U2AF1-S34','SRSF2-P95'],tandem['SRSF2-P95','U2AF1-S34']
print countgrp['U2AF1-S34'],countgrp['SRSF2-P95']


for j in grplist:
    export.write(j)
    for ij in grplist:
        val=0.0
        if j==ij:
            export.write('\t'+'0')
            continue
        if tandem[j,ij]==0:
            export.write('\t'+'0')
            continue
        if countgrp[j]>countgrp[ij]:
            val=float(tandem[j,ij])/float(countgrp[ij])
            export.write('\t'+str(val))
            continue
        else:
            val=float(tandem[j,ij])/float(countgrp[j])
            export.write('\t'+str(val))
            continue
            
    export.write('\n')
        

    