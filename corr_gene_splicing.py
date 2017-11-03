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

from multiprocessing import Manager,Process

lis=[]
tandem = dict()
dem=dict()
new=dict()
samplelis=[]
key_length=[]
list_g=[]
lis1=[]

lis2=[]

def worker(i):
  """thread worker function"""
  print 'Worker:'+str(i)
  ki=(i*j)
  q=0
  name='/Users/meenakshi/Documents/AMLComplete/LAML1/AltResults/Unbiased/junctions/4781/correlation4781all_details_'+str(i)+'.txt'
  export_ob1=open(name,'w')	
  if ((i+1)*j)<len(lines):
	
	q=((i+1)*j)
  else:
        ki=i*j
	q=len(lines)
  print ki,q
  for l in range(ki,q):
   if l== 0:
	continue
   else:
    t=string.split(lines[l],'\t')
    t=t[9:]
    for k in range(l+1,len(lines)):
        list1=[]
        list2=[]
        ind=[]
        p=string.split(lines[k],'\t')
        p=p[9:]
        for i in range(len(t)-1):
            if(t[i]!='' and p[i]!=''):
                ind.append(i)
            else:
                continue
        for i in range(len(ind)-1):
            list1.append(float(t[ind[i]]))
            list2.append(float(p[ind[i]]))
        
        if len(list1)==0 or len(list2)==0:
            #print l
	    correc[gene_label[l],gene_label[k]]=0
            continue
        else:
            max_value=max(list1)
            min_value=min(list1)
            me_std=np.std(list1)
	    
            r=sum(1 for i in list1 if i > max_value-me_std)
            s=sum(1 for i in list1 if i < min_value+me_std)
            max_value=max(list2)
            min_value=min(list2)
            me_std1=np.std(list2)
            r1=sum(1 for i in list2 if i > max_value-me_std1)
            s1=sum(1 for i in list2 if i < min_value+me_std1)
            #print r,s,s1,r1
            if (r > 2 and s > 2 and r1 >2 and s1 > 2):
                coefr=pearsonr(list1,list2)
                coef=coefr[0]
		if me_std==0 or me_std1==0:
		   correc[gene_label[l],gene_label[k]]=0
		else:
		    correc[gene_label[l],gene_label[k]]=coef
		
		pvalu=coefr[1]
                if ((coef > 0.15 or coef <-0.15)):
		    #print gene_label[l],gene_label[k]
		    if gene_label[l]==gene_label[k]:
			continue
		    else:
			export_ob1.write(gene_label[l])
			export_ob1.write('\t')
			export_ob1.write(gene_label[k])
			export_ob1.write('\t')
			export_ob1.write(str(coef))
			export_ob1.write('\n')
            else:
		    
		    coefr=pearsonr(list1,list2)
		    coef=coefr[0]
		    if me_std==0 or me_std1==0:
			correc[gene_label[l],gene_label[k]]=0
		    else:
			#correc[gene_label[l],gene_label[k]]=coef
			correc[gene_label[l],gene_label[k]]=0
  export_ob1.close()              
  return

def strip_first_col(fname, delimiter=None):
    with open(fname, 'r') as fin:
        for line in fin:
            
            try:
               yield line.split(delimiter, 1)[1]
            except IndexError:
               continue

def genelist(fname):
    for line in open(fname,'rU').xreadlines():
        line = line.rstrip(os.linesep)
        t=string.split(line,'\t')
	val=t[0]+' '+t[1]
        lis.append(val)
        #print t[0]
    return lis

def sample(fname):
    head=0
    for line in open(fname,'rU').xreadlines():
        line = line.rstrip(os.linesep)
	if head ==0:
	    t=string.split(line,'\t')
	    #print t
	    for p in range(9,len(t)):
		samplelis.append(t[p])
	    head=1
        else:
	    break;
    return samplelis

def create_corr_files():
    export_corrmat=open('/Users/meenakshi/Documents/AMLComplete/LAML1/AltResults/Unbiased/junctions/4781/correlation4781_matrixall.txt','w')
    #export_corrmat=open('correlation_matrix_up.txt','w')
    for i in range(5):
	read_genes=open('/Users/meenakshi/Documents/AMLComplete/LAML1/AltResults/Unbiased/junctions/4781/correlation4781all_details_'+str(i)+'.txt','r')
	line = read_genes.readlines()
        print line
        read_genes.close()
	for li in range(len(line)):
	    t=string.split(line[li],'\t')
	    temp=t[0]
            temp2=t[1]
            if temp not in tandem:
                        tandem[temp] = [temp2,]
            elif temp2 not in tandem[temp]:
              		tandem[temp].append(temp2)
            if temp2 not in tandem:
        		tandem[temp2] = [temp,]
            elif temp not in tandem[temp2]:
              		tandem[temp2].append(temp)
    export_ob.close() 
    for ke in tandem:
	export_corrmat.write('\t'+ke)
    export_corrmat.write('\n')
    for k1 in tandem:
	export_corrmat.write(k1+'\t')
	for k2 in tandem:
	  #if(k2 in tandem[k1]):
	    if k1==k2:
		export_corrmat.write(str(1)+'\t')
	    else:
		try:
		    export_corrmat.write(str(correc[k1,k2])+'\t')
		except KeyError:
		    export_corrmat.write(str(correc[k2,k1])+'\t')
	  #else:
	   # export_corrmat.write(str(0)+'\t')
	export_corrmat.write('\n')
	    #export_corrmat.write('\n')
    export_corrmat.close()
    invalid=dict()    
    invalid[1]=[]   
    count=0  
    for key in tandem:
    
        export_object.write(key+'\t'+str(tandem[key])+'\n')
        #print str(tandem[key])
        rem=[]
        if key in invalid[1]:
            continue
        else:

             invalid[1].append(key)
             dem[count]=[key,]
             
             for val in tandem[key]:
                 if val not in dem[count]:
                     dem[count].append(val)
                     if len(tandem[val]) >1:
                         for k in tandem[val]:
			        
	                         if k not in dem[count] and key not in k:
                                     dem[count].append(k)
				 else:
				    invalid[1].append(k)
				 
		     else:
                         if tandem[val] not in dem[count] and key not in tandem[val]:
                               dem[count].append(tandem[val])
			 else:
			    invalid[1].append(tandem[val])
	values=dem.values()
	dem[count]=sorted(dem[count])
	if dem[count] in values:
	     dem[count]=[]
	if len(dem[count])>2:		
	    count=count+1
	else:
	    dem[count]=[]
	
	
    for key in dem:
             print str(key)
             export_object2.write(str(key)+'\t'+str(dem[key]))
             export_object2.write('\t'+str(len(dem[key])))
             export_object2.write('\n')
	     #key_length[int(key)]=len(dem[key])
    export_object.close()
    export_object2.close()
	    

def create_sub_file():
    export_ob3=open("/Users/meenakshi/Documents/AMLComplete/LAML1/AltResults/Unbiased/junctions/4781/alt_spliced4781.txt",'w')
    head=0
    for line in open('/Users/meenakshi/Documents/AMLComplete/LAML1/AltResults/Unbiased/junctions/4781/top_Alt_sig4781.txt','rU').xreadlines():
        if head==0:
            
            k=string.split(line,'\t')
            count1 = len(k)-9
            #count1=int(count1/3)
            count1=int(count1)
            export_ob3.write(line)
            #print line
	    #export_ob3.write('\n')
            head=1
        else:
            #print count1
            t=string.split(line,'\t')
            if(len(t)-9==count1):
                #print line
                export_ob3.write(line)
             #   export_ob3.write('\n')
    export_ob3.close()
    #sys.exit()
    
   
	
def strip_first_col(fname, delimiter=None):
    with open(fname, 'r') as fin:
        for line in fin:
            try:
               yield line.split(delimiter, 1)[1]
            except IndexError:
               continue


def splitfile():
    count =0
    head=0
    export_object = open("/Volumes/MyPassport/Figure1_final/exp.Leucecorr_1.txt",'w')
   

    for line in open('/Volumes/MyPassport/Figure1_final/exp.Leucecorr.txt','rU').xreadlines():
        if head==0:
         
            head=1
            continue
       
        line=line.rstrip(os.linesep)
        line1=string.split(line,'\t')
        line1=line1[1:]
        for i in line1:
            if i=='0.0':
                continue
            if abs(float(i))>0.4:
                export_object.write(line)
                export_object.write('\n')
                count=count+1
    print count
    
if __name__ == '__main__':
    #manager=Manager()
    #start = timeit.default_timer()
    #correc=manager.dict()
   # create_sub_file()
   
    splitfile();sys.exit()   
    geneexp= defaultdict(list)
    genelst=[]
    totalgene=[]
    splicelst=[]
    test_name="/Volumes/MyPassport/Figure1_final/exp.Leucesymb.txt"
    text_file = open(test_name,'r')
    lines = text_file.readlines()
    text_file.close()
    export_object2 = open("/Volumes/MyPassport/Figure1_final/exp.Leucecorr.txt",'w')
    coef={}
    #create_sub_file()
    head=0
    for i in lines:
        if head==0:
            head=1
            continue
        i = i.rstrip(os.linesep)
        i=string.split(i,'\t')
        valu=i[1:]
        valu=[float(ij) for ij in valu]
       
        geneexp[str(i[0])]=valu
        genelst.append(str(i[0]))
      
        
    head=0
    for line in open('/Volumes/MyPassport/Figure1_final/exp.splicing.txt','rU').xreadlines():
        if head==0:
            head=1
            continue
        line=line.rstrip(os.linesep)
        line=string.split(line,'\t')
        gene=string.split(line[0],'&')[0]
        
        val=line[1:]
        val=[float(i) for i in val]
        
        try:
           if gene not in totalgene:
            totalgene.append(gene)
           coefr=pearsonr(val,geneexp[gene])
          
           coef[line[0],gene]=round(coefr[0],3)
           
           
        except Exception:
            coef[line[0],gene]=0.0
        
        splicelst.append(line[0])
        
    print len(totalgene)
        
    for j in totalgene:
       
      
        export_object2.write('\t'+j)
    export_object2.write('\n')
        
    for i in splicelst:
        export_object2.write(i)
        for j in totalgene:
            try:
                export_object2.write('\t'+str(coef[i,j]))
            except KeyError:
                export_object2.write('\t'+'0.0')
        export_object2.write('\n')
        
    
    
#T=np.loadtxt("C:/Users/venz6v/Documents/classes/Research/Cardiac_SingleCell_Seq/trial.txt")k]
    
            
#    export_ob = open("/Users/meenakshi/Documents/AMLComplete/LAML1/AltResults/Unbiased/junctions/4781/corr_results_splicing_fl4781.txt",'w')
#    export_object = open("/Users/meenakshi/Documents/AMLComplete/LAML1/AltResults/Unbiased/junctions/4781/corr_res4781_fl.txt",'w')
#    export_object2 = open("/Users/meenakshi/Documents/AMLComplete/LAML1/AltResults/Unbiased/junctions/4781/corr_res14781_fl.txt",'w')
#    #export_in=open("/Users/meenakshi/Documents/AMLComplete/LAML1/AltResults/Unbiased/junctions/4781/events.gff",'w')
#    gene_label=genelist(test_name)
#    sample_list=sample(test_name)
#    
#    count=0
#    jobs = []
#    j=int(len(lines)/5)
#   
##    processes = [multiprocessing.Process(target=worker, args=(x,)) for x in range(5)]
##    for p in processes:
##	p.start()
##    for p in processes:
##	p.join()
#
#    create_corr_files()
#  
#    stop= timeit.default_timer()
#    print stop- start 

