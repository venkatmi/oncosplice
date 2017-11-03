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
dem=dict()
new=dict()
samplelis=[]
key_length=[]
list_g=[]
lis1=[]
correc=dict()
lis2=[]
gene_labels=[]
from multiprocessing import Manager,Process


def worker(i):
 """thread worker function"""
 print 'Worker:'+str(i)
 print len(lines),len(tfs)
 counter=0
 for l in range(len(lines)):
    k_keys_sorted_by_values={}
    #t=lines[l].rstrip(os.linesep)
    #t=string.split(t,'\t')
    t=lines[l]
    t1=gene_labels[l]
    #print t1
    t2=[]
    if t1 in tfs:
      export1=open(t1+str(i)+"_permute2.txt","w")
      
      for i in t:
         t2.append(int(i))
      ones=np.count_nonzero(t2)
      #print ones
      counter=counter+1
      pvaluelist=[]
      for iter in range(10):
         counter=counter+1
         newones=[]
         newarray=[]
         indices=[]
         correc1=dict()
         corind=dict()
         
         newones=np.random.choice(len(t),ones, replace=False)
         for n in range(len(t)):
            if n in newones:
               newarray.append(int(1))
            else:
               newarray.append(int(0))
         
         for k in range(len(lines)):
               #if k ==0:
                # continue
               indices=[]
               #list1=[]
               #list2=[]
               ind=[]
               p=lines[k]
            
               
               #for i in range(len(t)-1):
                #if(newarray[i]!='' and p[i]!=''):
                 #   ind.append(i)
                #else:
                 #   continue
               for i in range(len(p)-1):
                #list1.append(float(newarray[ind[i]]))
                #list2.append(float(p[ind[i]]))
                #if i in newones:
                if float(newarray[i])==float(1):
                  indices.append(float(p[i]))
               #if len(list1)==0 or len(list2)==0:
                #correc1[gene_labels[k-1]]=0
                #continue
               #else:
               coefr=pearsonr(p,newarray)
               coef=coefr[0]
	       correc1[gene_labels[k]]=coef
               corind[gene_labels[k]]=np.mean(indices)
               # print indices
         
         print counter 
         k_keys_sorted_by_values = heapq.nlargest(60,correc1, key=correc1.get)
       
         counter2=0
         ave=float(0)
         average=float(0)
         for key in k_keys_sorted_by_values:
           # print corind[key]
            #for lis in corind[key]:
                  ave=float(ave)+float(corind[key])
                  counter2=counter2+1
        # print str(ave),str(counter2)
         average=float(ave)/float(counter2)
      #  print str(average)
         print str(average)
         pvaluelist.append(float(average))
         
      for i in range(len(pvaluelist)):
       print (pvaluelist[i])
       export1.write(str(pvaluelist[i])+'\n')        
                 
 

def strip_first_col(fname, delimiter=None):
    with open(fname, 'r') as fin:
        for line in fin:
            
            try:
               yield line.split(delimiter, 1)[1]
            except IndexError:
               continue

def genelist(fname):
    head=0
    lis=[]
    array=[]
    ke=[]
    for line in open(fname,'rU').xreadlines():
        line = line.rstrip(os.linesep)
        t=string.split(line,'\t')
         
	val=t[0]
        if head==0:
         
         if t[0]=='AltAnalyze_ID':
             continue
         else:
           lis.append(val)
         head=1
         continue
        
        ke=t[1:]
        newarray=[]
        if len(t)>1:
         for i in range(len(ke)):
          if ke[i]=="":
            newarray.append(float(0))
          else:
            newarray.append(float(str(ke[i])))
         
          #print ke[i]
         array.append(newarray)
         lis.append(val)
        
    #print lis
 
    return lis,array

def sample(fname):
    head=0
    for line in open(fname,'rU').xreadlines():
        line = line.rstrip(os.linesep)
	if head ==0:
	    t=string.split(line,'\t')
	    #print t
	    for p in range(1,len(t)):
		samplelis.append(t[p])
	    head=1
        else:
	    break;
    return samplelis

def create_corr_files(correc,filename, tfs, gene_labels):
    export_corrmat=open(filename[:-4]+'-corr.txt','w')
    #export_corrmat=open('correlation_matrix_up.txt','w')
    for li in range(len(tfs)):
        export_corrmat.write('\t'+tfs[li])
    export_corrmat.write('\n')
    for i in range(len(gene_labels)):
        export_corrmat.write(gene_labels[i]+'\t')
        
	for li in range(len(tfs)):
            try:
	        export_corrmat.write(str(correc[gene_labels[i],tfs[li]])+'\t')
	    except Exception:
               print traceback.format_exc()
               export_corrmat.write('NA')
	  #else:
	   # export_corrmat.write(str(0)+'\t')
	export_corrmat.write('\n')
	    #export_corrmat.write('\n')
    export_corrmat.close()
   
def strip_first_col(fname, delimiter=None):
    with open(fname, 'r') as fin:
        for line in fin:
            try:
               yield line.split(delimiter, 1)[1]
            except IndexError:
               continue

def runTFCorrelationAnalysis(query_exp_file,query_tf_file):
    query_data = open(query_exp_file,'rU')
    lines = query_data.readlines()
    print "Number or rows in file:",len(lines)
    query_data.close()

    tfs=genelist(query_tf_file)
    gene_labels=genelist(query_exp_file)
    
    correc=create_corr_matrix(lines,tfs,gene_labels)
    #create_corr_files(correc,query_exp_file,tfs,gene_labels)

if __name__ == '__main__':
    import getopt
    filter_rows=False
    filter_file=None
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
      print "Insufficient arguments";sys.exit()
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','t='])
        #print sys.argv[1:]
        for opt, arg in options:
            if opt == '--i': query_exp_file=arg
            elif opt == '--t': query_tf_file=arg
    #query_data = open(query_exp_file,'rU')
    #lines = query_data.readlines()
    #print "Number or rows in file:",len(lines)
    #query_data.close()
    
    tfs,lines=genelist(query_tf_file)
    print tfs
    gene_labels,lines=genelist(query_exp_file)
    
    print "Number or rows in file:",len(lines)
    processes = [multiprocessing.Process(target=worker, args=(x,)) for x in range(1)]
    for p in processes:
	p.start()
    for p in processes:
	p.join()
   # correc=create_corr_matrix(lines,tfs,gene_labels) 
   # runTFCorrelationAnalysis(query_exp_file,query_tf_file)