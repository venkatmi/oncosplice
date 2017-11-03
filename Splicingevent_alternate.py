#!/usr/bin/env python

import numpy as np
import sys,string
import os
import os.path
from numpy import corrcoef, sum, log, arange
from scipy.stats.stats import pearsonr
import traceback
import heapq
import multiprocessing
tandem = dict()
lis=[]



fname="/Users/meenakshi/Documents/April12_updates/Alt_events_all.txt"
for line in open(fname,'rU').xreadlines():
        
    line = line.rstrip(os.linesep)
    t=string.split(line,'\t')
    val=t[0]
    lis.append(val)
    tandem[val]=t[1]
    #print lis
 
fname="/Users/meenakshi/Documents/April12_updates/exp.U2AF1like.txt"

export=open("/Users/meenakshi/Documents/PSIAnnotFile_commonevents2.txt","a")
head=0
count=0


for line in open(fname,'rU').xreadlines():
    line = line.rstrip(os.linesep)
    gene_val=[]
    cor=[]
    if head==0:
	head=1
        continue
    else:
        t=string.split(line,'\t')
        if t[3] in lis:
            for i in t[4:len(t)]:
                gene_val.append(float(i))
            head2=0
            for lin in open(fname,'rU').xreadlines():
                gene_val1=[]
                lin = lin.rstrip(os.linesep)
                if head2==0:
                    head2=1
                else:
                    t1=string.split(lin,'\t')
                    for ii in t1[4:len(t1)]:
                        gene_val1.append(float(ii))
                        
                    co=np.corrcoef(gene_val,gene_val1)[0][1]
                    cor.append(co)
                    if t1[3]==tandem[t[3]]:
                        co1=co
                        
                            
            cor = sorted(cor, key=float)
            k=cor.index(co1)
            count+=1
            #print t[3],str(co1),str(cor[0]),str(k)
            print count
            export.write(t[3]+'\t'+str(k)+'\n')
            
            
            
                        
            
                            
                            
                        
                        
                        
                    
                
       
        
            
        