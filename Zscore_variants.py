#!/usr/bin/env python

import os
import os.path
global test_list
import matplotlib.pyplot as plt
import copy
from collections import Counter
from collections import defaultdict
import string
tandem = {}
head=0
import math
def zscore(r,n,N,R):               
    r = float(r)       #number of patients with the gene in that cluster
    n = float(n)       #number of patients in that cluster
    N = float(N)                     #total number of patients
    R = float(R)      #number of patients that have that gene in dataset
    if (R-N) == 0: return 0
    elif r==0 and n == 0: return 0
    else:
        try:
            #try:
            z = (r - n*(R/N))/math.sqrt(n*(R/N)*(1-(R/N))*(1-((n-1)/(N-1))))
            print z
            return z
            #except ZeroDivisionError: print 'r,n,R,N: ', r,n,R,N;kill
        except Exception:return 0

def initiaterec():
    text_file = open("/Volumes/MyPassport/Leucegene/combined_all_samples/High_inputfile-output.txt", "rU")
    
    line = text_file.readlines()

    text_file.close()
    export=open('/Volumes/MyPassport/Leucegene/combined_all_samples/High_inputfile-output_zscore.txt','w')
    head=0
    a = defaultdict(list)
    export.write('UID'+'\t'+'PRPF40B-2'+'\t'+'No.of terms in group'+'\t'+'PML-RARA'+'\t'+'No.of terms in group'+'\t'+'SF3B1'+'\t'+'No.of terms in group'+'\t'+'U2AF1-like'+'\t'+'No.of terms in group'+'\t'+'ZRSR2'+'\t'+'No.of terms in group'+'\t'+'CBFB-MYH11'+'\t'+'No.of terms in group'+'\t'+'SRSF2'+'\t'+'No.of terms in group'+'\t'+'U2AF2'+'\t'+'No.of terms in group'+'\t'+'SRSF2-like'+'\t'+'No.of terms in group'+'\t'+'PRPF40B'+'\t'+'No.of terms in group'+'\t'+'HNRNPK'+'\t'+'No.of terms in group'+'\t'+'SF3A2'+'\t'+'No.of terms in group'+'\t'+'RBM38'+'\t'+'No.of terms in group'+'\t'+'THRAP3'+'\t'+'No.of terms in group'+'\t'+'MLL'+'\t'+'No.of terms in group'+'\t'+'others'+'\t'+'No.of terms in group'+'\t'+'NPM1'+'\t'+'No.of terms in group'+'\t'+'SNRNP40'+'\t'+'No.of terms in group'+'\t'+'SRSF2-low'+'\t'+'No.of terms in group'+'\t'+'others-cl4'+'\t'+'No.of terms in group'+'\t'+'others-cl3'+'\t'+'No.of terms in group'+'\t'+'others-cl2'+'\t'+'No.of terms in group'+'\t'+'others-cl1'+'\t'+'No.of terms in group'+'\t'+'PRPF40B-2'+'\t'+'No.of terms in group'+'\t'+'RUNX1'+'\t'+'No.of terms in group'+'\t'+'U2AF1-A'+'\t'+'No.of terms in group'+'\t'+'U2AF1-B'+'\t'+'No.of terms in group'+'\t'+'total'+'\n')
    header=0
    for l in line:
        R=0
        l=l.rstrip('\r\n')
     
        t = string.split(l,'\t')
        N=len(t)-1
        
        n = defaultdict(int)
        r = defaultdict(int)
        #print a
        for i in range(1,len(t)):
            if head==0:
                grpname=string.split(t[i],":")[1]
                a[grpname].append(i)
                
                
        

            else:
                if float(t[i])>0.0:
                #if t[i]!="./.":
                    R+=1
                    for key in a:
                        if i in a[key]:
                            r[key]=r[key]+1
            for key in a:
                n[key]=len(a[key])
        
        if head==0:
            head=1
            continue
        export.write(t[0]+'\t')
        for key in a:
           # print N,R
            print t[0],str(r[key]),str(n[key]),N,R
            z=zscore(r[key],n[key],N,R)
                
            export.write(str(z)+'\t'+str(r[key])+'\t')
        export.write(str(R)+'\n')
    if header==0:
        for key in a:
            print key

initiaterec()
              
      
