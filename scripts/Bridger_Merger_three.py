#!/usr/bin/env python

#!/usr/bin/env python
import numpy as np
import sys,string
import os
import os.path
from collections import defaultdict
import sampleIndexSelection
import export
from bs4 import BeautifulSoup
import urllib2
import csv
import sys
import re
import math

Newlist=defaultdict(list)
Newval={}
genelst=defaultdict(list)
Allcomp={}

MergedOutput="/Volumes/Pass/MotifAnalyses/Bridger/Exons_MotifAnalyses/merged_output_allpvalues_nofold.txt"
output="/Volumes/Pass/MotifAnalyses/Bridger/Exons_MotifAnalyses/merged_output_allpvalues_nofold_upd.txt"
output1=open(output,"w")

for lin in open(MergedOutput,'rU').xreadlines():
    genes=[]
    s=lin.rstrip('\r\n')
    
    s1=string.split(s,'\t')
    sig=s1[0]
    
    if s1[2]=="GE":
        genes=[s1[1]]
        
        
    else:
        genes=string.split(s1[1],":")
    
    tool=s1[2]
    
    if 'Cisbp_denovo' in tool:
        tool="Cisbp_denovo"
    if "UpstreamIntron_known" in sig:
        sig = string.replace(sig,"UpstreamIntron_known","Upstream")
      
    if "Intron_known" in s1[0]:
        sig=string.replace(sig,"Intron_known","Combined_intron_new")
        
    if "Exons_known" in s1[0]:
        sig=string.replace(sig,"Exons_known","Exon")

    if "DownstreamIntron_known" in s1[0]:
        sig=string.replace(sig,"DownstreamIntron_known","Downstream")

        
    for i in range(len(genes)):
        if tool not in genelst[sig,genes[i].upper()]:
            genelst[sig,genes[i].upper()].append(tool)
            
        Newval[sig,tool,genes[i].upper()]=float(s1[3])
        if tool=="GE":
            sig1="Exon:"+sig
            Newval[sig1,tool,genes[i].upper()]=float(s1[3])
            
            genelst[sig1,genes[i].upper()].append(tool)
            sig1="Combined_intron_new:"+sig
            Newval[sig1,tool,genes[i].upper()]=float(s1[3])
            genelst[sig1,genes[i].upper()].append(tool)
        

for sig,genes in genelst:
    tools=[]
    cisbpact=True
    cisbpden=True
    tools=genelst[sig,genes]
   # if genes=="MBNL1":
   # print tools,sig
    a=len(tools)
    if 'Cisbp_Actual' in tools and 'Cisbp_denovo' in tools:
        a=a-1
        if Newval[sig,"Cisbp_Actual",genes]<Newval[sig,"Cisbp_denovo",genes]:
            cisbpden=False
        else:
            cisbpact=False
            
    pval=0.0
    count=0
    if a>1:
        if "Cisbp_Actual" in tools and cisbpact:
            output1.write(sig+"\t"+genes+"\t"+"Cisbp_Actual"+"\t"+str(Newval[sig,"Cisbp_Actual",genes])+"\t")
            count+=1
           # print str(Newval[sig,"Cisbp_Actual",genes])
            pval=pval-math.log10(Newval[sig,"Cisbp_Actual",genes])
        else:
            output1.write(sig+"\t"+genes+"\t"+"Cisbp_Actual"+"\t"+"NA"+"\t")
        if 'Cisbp_denovo' in tools and cisbpden:
            count+=1
            #print str(Newval[sig,"Cisbp_denovo",genes])
            pval=pval-math.log10(Newval[sig,"Cisbp_denovo",genes])
            output1.write(sig+"\t"+genes+"\t"+"Cisbp_denovo"+"\t"+str(Newval[sig,"Cisbp_denovo",genes])+"\t")
        else:
            output1.write(sig+"\t"+genes+"\t"+"Cisbp_denovo"+"\t"+"NA"+"\t")
            
        if "Clipseq" in tools:
            count+=1
            #print str(Newval[sig,"Clipseq",genes])
            pval=pval-math.log10(Newval[sig,"Clipseq",genes])
            output1.write(sig+"\t"+genes+"\t"+"Clipseq"+"\t"+str(Newval[sig,"Clipseq",genes])+"\t")
        else:
           output1.write(sig+"\t"+genes+"\t"+"Clipseq"+"\t"+"NA"+"\t")
        if "GE" in tools:
            count+=1
            #print str(Newval[sig,"GE",genes])
            pval=pval-math.log10(Newval[sig,"GE",genes])
            output1.write(sig+"\t"+genes+"\t"+"GE"+"\t"+str(Newval[sig,"GE",genes])+"\t")
        else:
            output1.write(sig+"\t"+genes+"\t"+"GE"+"\t"+"NA"+"\t")
        pval=pval/float(count)
        output1.write(str(pval)+"\n")
            
        
        
            
            