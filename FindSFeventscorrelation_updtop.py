#!/usr/bin/env python
import os
import os.path
global test_list
#import matplotlib.pyplot as plt
#import copy
from collections import Counter
from collections import defaultdict
import string
import numpy
import scipy
from scipy import stats
tandem = {}
head=0
import math


def findcorrelations():
    # input psi file  
    SplicingeventsFile="round2_input-filtered.txt"
    # filtered expression file for selected genes
    SFexp_file = open("Tophat_Exp.txt", "rU")
    lines = SFexp_file.readlines()
    SFexp_file.close()
    
    #writing the splicing factor results
    export=open('Tophat_SF_comps.txt','w')
    
    #processing the splicing events file
    head=True
    counter=0
    for exp in lines:
        exp=exp.rstrip('\r\n')
        if head:
            export.write('Splicing Factor'+'\t'+'events count'+'\n')
            
            head = False
            continue
        firstLine=True
        s = string.split(exp,'\t')
        Corrsflist=[]
        sfname=s[0]
        count=0
        
       # export=open('/Users/meenakshi/Documents/ICGS_guidefiles/Guide3_firstroundicgs-correlations.txt','w')
        sfvalues = map(lambda x: float(x), s[1:])
        for splevent in open(SplicingeventsFile,'rU').xreadlines():
            splevent=splevent.rstrip('\r\n')
            list1=[]
            list2=[]
            eve = string.split(splevent,'\t')
            if firstLine:
               firstLine = False
               continue
            else:
               
                    
                    for val in range(1,len(eve)):
                        
                        try:
                            eve[val]=float(eve[val])
                            list1.append(float(eve[val]))
                            list2.append(float(sfvalues[val-1]))
                        except Exception:
                            continue
                            
                        
                        
               # values = numpy.ma.masked_values(values,0.000101)
            coefr=scipy.stats.pearsonr(list1,list2)[0]
           
            if abs(coefr)>0.5:
                    count+=1
                    Corrsflist.append([eve[0],coefr])
        
        export.write(sfname+'\t'+str(count)+'\n')
        filename=SplicingeventsFile[:-4]+"_"+sfname+"_"+str(count)+"tophat.txt"
        if count>20:
            exportgene=open(filename,"w")
            exportgene.write("SplicingEvent"+"\t"+"Code"+"\n")
        
            for si in range(len(Corrsflist)):
                exportgene.write(Corrsflist[si][0]+"\t"+"Ae"+"\t"+str(Corrsflist[si][1])+"\n")
        
        counter+=1
        print counter

findcorrelations()
