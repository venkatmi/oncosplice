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


fname="/Users/meenakshi/Documents/April12_updates/MergedFiles.txt"
fexp="/Users/meenakshi/Documents/April12_updates/Leuce_exp_log-Symbol-filtered.txt"
export=open("/Users/meenakshi/Documents/exp.splicing-pairwise_correlations.txt","w")
export1=open("/Users/meenakshi/Documents/exp.splicing-expression.txt","w")
head=0

prev_gene=""
gene_val=[]
val=[]
mat=[]
counter=1

val_mean=[]
me=0.0
for line in open(fname,'rU').xreadlines():
    line = line.rstrip(os.linesep)
    found=0
    if head==0:
	head=1
    else:
        counter+=1
        t=string.split(line,'\t')
        for i in t[4:len(t)]:
            
            val.append(float(i))
        #if t[0]=="BANK1":
                
        if t[0]==prev_gene or prev_gene=="":

            gene_val.append(val)
            
        else:
            exp_cor=[]
            #val_mean=np.median(gene_val,axis=0)
            for lin in open(fexp,"rU").xreadlines():
                lin = lin.rstrip(os.linesep)
                t1=string.split(lin,'\t')
                if prev_gene==t1[0]:
                    val1=[]
                    for i in t1[1:len(t1)]:
                        val1.append(float(i))
                    for vl in gene_val:
                        try:
                            exp_cor.append(np.corrcoef(vl,val1)[0][1])
                            found=1
                        except:
                            exp_cor.append(0.0)
                            found=1
                    break
            if found==1:
                export1.write(prev_gene)
                for ex in exp_cor:
                    export1.write('\t'+str(ex))
                export1.write('\n')
          
                
                    
          
            try:
                me_list=[]
                me=float(np.mean(list(np.corrcoef(gene_val)[np.triu_indices_from(np.corrcoef(gene_val),k=1)])))
                me_list=list(np.corrcoef(gene_val)[np.triu_indices_from(np.corrcoef(gene_val),k=1)])
                export.write(prev_gene)
                for ls in me_list:
                    export.write('\t'+str(ls))
                export.write('\n')
                #med=float(np.median(list(np.corrcoef(gene_val)[np.triu_indices_from(np.corrcoef(gene_val),k=1)])))
                #absmean=float(np.mean(list(abs(np.corrcoef(gene_val)[np.triu_indices_from(np.corrcoef(gene_val),k=1)]))))
                #absmed=float(np.median(list(abs(np.corrcoef(gene_val)[np.triu_indices_from(np.corrcoef(gene_val),k=1)]))))
            #print str(me),str(med),str(absmean);sys.exit()
                
            except Exception:
                continue
                me=float(0)
                med=float(0)
                absmean=float(0)
                absmed=float(0)
            #tandem[prev_gene]=[me,med,absmean,absmed,exp_cor]
            gene_val=[]
            gene_val.append(val)
       
        prev_gene=t[0]
        
        val=[]
        print counter
        
#for gene in tandem:
 #   export.write(gene+'\t'+str(tandem[gene][0])+'\t'+str(tandem[gene][1])+'\t'+str(tandem[gene][2])+'\t'+str(tandem[gene][3])+'\t'+str(tandem[gene][4])+'\n')
        