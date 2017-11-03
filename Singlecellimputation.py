#!/usr/bin/env python

import os
import os.path
global test_list
import matplotlib.pyplot as plt
import copy
from collections import Counter
from collections import defaultdict
import string
import numpy as np
tandem = {}
import sys
head=0
import math


def initiaterec():
    genelist=[]
    grpsdict=defaultdict(list)
    text_file1 = open("/Users/meenakshi/Downloads/hscardiacrnaseq/ExpressionInput/filteredExp.PSI_file-filtered.txt","rU")
    lines = text_file1.readlines()
    text_file1.close()
    
    export=open('/Users/meenakshi/Downloads/hscardiacrnaseq/ExpressionInput/filteredExp.PSI_file-filtered_imputed.txt','w')
    head=0
    for line in lines:
            lin=line.rstrip('\r\n')
            lin=string.split(lin,"\t")
            if head==0:
                export.write(line)
                for i in range(len(lin)):
                    if ':' in lin[i]:
                        grps=string.split(lin[i],":")
                        grpsdict[grps[0]].append(i)
                head=1
                continue
            newlin=[]
            
            grp_means={}
            a=[]
            export.write(lin[0]+"\t")
            
            for grp in grpsdict:
                a=[lin[x] for x in grpsdict[grp]]
                values=[]
                for value in a:
                    try: values.append(float(value)) 
                    except Exception: values.append(0.000101)
                
                values = np.ma.masked_values(values,0.000101)
                grp_means[grp]=np.ma.mean(values)
           
                
            for i in range(1,len(lin)):
                if lin[i]==0 or lin[i]=="":
                    for grp in grpsdict:
                        if i in grpsdict[grp]:
                            if grp_means[grp] =="__":
                                imputemean="0.0"
                            else:
                                imputemean=grp_means[grp]
                        
                    newlin.append(str(imputemean))
                else:
                    
                    newlin.append(lin[i])
                
            
            export.write("\t".join(newlin))
            export.write("\n")
initiaterec()
