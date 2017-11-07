#!/usr/bin/env python
import numpy as np
import sys,string
import os
import os.path
from numpy import corrcoef, sum, log, arange
from scipy.stats.stats import pearsonr
from collections import Counter
from collections import defaultdict
import traceback
tandem = defaultdict(int)
dem=dict()
lst=[]
delt=[]
strlst=[]

spliceconn="/Volumes/MyPassport/Leucegene_U2AF1vsU2AF1like/correlation_0.8_upd.txt"
export1=open("/Volumes/MyPassport/Leucegene_U2AF1vsU2AF1like/correlation_0.8_upd_nodup.txt",'w')

for line in open(spliceconn,'rU').xreadlines():
    line1 = line.rstrip(os.linesep)
  #  line=line.replace('\"','')
    t=string.split(line1,'\t')
    
    gene1=t[0]
    gene2=t[1]
    st=gene1+gene2
    if st not in strlst:
        strlst.append(st)
        export1.write(line)
strlst=[]
    

spliceconn1="/Volumes/MyPassport/Leucegene_U2AF1vsU2AF1like/correlation_0.8_upd_nodup.txt"
component="/Volumes/MyPassport/Leucegene_U2AF1vsU2AF1like/Results_components/Components_0.8_upd.txt"

cnt=Counter()




for line in open(spliceconn1,'rU').xreadlines():
    line = line.rstrip(os.linesep)
  #  line=line.replace('\"','')
    t=string.split(line,'\t')
    
    gene1=t[0]
    gene2=t[1]
    
    #print gene,val;sys.exit()
    if gene1!=gene2 and gene2 not in delt:

        #print gene1,gene2
        
        
        cnt[gene1]+=1
        if gene1 not in delt:
            delt.append(gene1)



print cnt["54761"]

export=open("density_0.8_upd2.txt",'w')

for line1 in open(component,'rU').xreadlines():
  
    line1 = line1.rstrip(os.linesep)
   # line1=line1.replace('\"','')
    t=string.split(line1,'\t')
    event=t[0]
    comp=t[1]
  #  t[2]=t[2].replace("-",".")
   # t[2]=t[2].replace("&",".")
    #splice2=t[2]
    try:
    #print splice1,splice2;sys.exit()
        #if event=='1':
         #   print cnt[event];sys.exit()
        tandem[comp]=tandem[comp]+cnt[event]
        if event=='6':
            print tandem[comp]
    except Exception:
        continue


for val in tandem:
    export.write(str(val)+'\t'+str(tandem[val])+'\n')

