#!/usr/bin/env python

#!/usr/bin/env python
import numpy as np
import pylab as pl
import sys,string
import os
import os.path
from collections import defaultdict
from sklearn.cluster import KMeans

import statistics
import random
import UI
import export; reload(export)
import re
import fishers_exact_test
import traceback
import warnings
import math

def FishersExactTest(r,n,R,N):
    z=0.0
    """
    N is the total number of genes measured (Ensembl linked from denom) (total number of ) (number of exons evaluated)
    R is the total number of genes meeting the criterion (Ensembl linked from input) (number of exonic/intronic regions overlaping with any CLIP peeks)
    n is the total number of genes in this specific MAPP (Ensembl denom in MAPP) (number of exonic/intronic regions associated with the SF)
    r is the number of genes meeting the criterion in this MAPP (Ensembl input in MAPP) (number of exonic/intronic regions with peeks overlapping with the SF)
    
    With these values, we must create a 2x2 contingency table for a Fisher's Exact Test
    that reports:
     
    +---+---+    a is the # of IDs in the term regulated
    | a | b |    b is the # of IDs in the term not-regulated 
    +---+---+    c is the # of IDs not-in-term and regulated
    | c | d |    d is the # of IDs not-in-term and not-regulated
    +---+---+
    
    If we know r=20, R=80, n=437 and N=14480
    +----+-----+    
    | 20 | 417 |  437   
    +----+-----+    
    | 65 |13978| 14043
    +----+-----+
      85  14395  14480
    """
    if (R-N) == 0: return 0
    elif r==0 and n == 0: return 0
    else:
        try:
            #try:
            z = (r - n*(R/N))/math.sqrt(n*(R/N)*(1-(R/N))*(1-((n-1)/(N-1))))
          
            #except ZeroDivisionError: print 'r,_n,R,N: ', r,_n,R,N;kill
        except Exception: print (r - n*(R/N)), n*(R/N),(1-(R/N)),(1-((n-1)/(N-1))),r,n,N,R;kill
    a = r; b = n-r; c=R-r; d=N-R-b
    table = [[int(a),int(b)], [int(c),int(d)]]

    """
    print a,b; print c,d
    import fishers_exact_test; table = [[a,b], [c,d]]
    ft = fishers_exact_test.FishersExactTest(table)
    print ft.probability_of_table(table); print ft.two_tail_p()
    print ft.right_tail_p(); print ft.left_tail_p()
    """
    
    try: ### Scipy version - cuts down rutime by ~1/3rd the time
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",category=RuntimeWarning) ### hides import warnings
            oddsratio, pvalue = stats.fisher_exact(table)
       # print pvalue
        return pvalue,z
    except Exception:
        ft = fishers_exact_test.FishersExactTest(table)
       # print ft.two_tail_p()
        return ft.two_tail_p(),z
    

            
def header_file(fname, delimiter=None):
    head=0
    header=[]
    newheader=[]
    with open(fname, 'rU') as fin:
        for line in fin:
            line = line.rstrip(os.linesep)
            header=string.split(line,'\t')
        
            if len(header)>2:
                if head==0:
                    
                    for i in range(1,len(header)):
                        iq=string.split(header[i],".")[0]
                        newheader.append(iq)
                    head=1
                else:break
            else:
               
                if header[0] not in newheader:
                    newheader.append(header[0]) 
    print len(newheader)           
    return newheader

def Enrichment(Guidefile,mutdict,mutfile,Expand,header):
    
    X=defaultdict(list)
    prev=""
    head=0
    group=defaultdict(list)
    mut=export.findFilename(mutfile)
    exportnam=Guidefile[:-4]+mut[:-4]+'enrichment.txt'
    export_enrich=open(exportnam,"w")
    export_enrich.write("Mutations"+"\t"+"Cluster"+"\t"+"Pvalue"+"\t"+"r"+"\t"+"R"+"\t"+"n"+"\t"+"z-score"+"\t"+"Sensitivity"+"\t"+"Specificity"+"\n")
    if Expand=="yes":
        header2=header_file(Guidefile)
        
        for line in open(Guidefile,'rU').xreadlines():
            if head >0:
                line=line.rstrip('\r\n')
                q= string.split(line,'\t')
                for i in range(1,len(q)):
                    if q[i]==str(1):
                        group[q[0]].append(header2[i-1])
           
            else:
                head+=1
                continue
    else:
        for line in open(Guidefile,'rU').xreadlines():
            line=line.rstrip('\r\n')
            line=string.split(line,'\t')
            #for i in range(1,len(line)):
            group[line[2]].append(line[0])
      
    for kiy in mutdict:
        
        groupdict={}
        remaining=[]
        remaining=list(set(header) - set(mutdict[kiy]))
        groupdict[1]=mutdict[kiy]
        groupdict[2]=remaining
        for key2 in group:
           
               
            r=float(len(group[key2])-len(list(set(group[key2]) - set(mutdict[kiy]))))
            n=float(len(group[key2]))
            R=float(len(set(mutdict[kiy])))
            N=float(len(header))
            #print kiy,key2,r,n,R,N
            if r==0:
                pval=float(1)
                z=float(0)
            else:
                try:
                    pval,z=FishersExactTest(r,n,R,N)
                    export_enrich.write(str(kiy)+"\t"+str(key2)+"\t"+str(pval)+"\t"+str(r)+"\t"+str(R)+"\t"+str(n)+"\t"+str(z)+"\t"+str(float(r)/(float(R)))+"\t"+str(float(r)/(float(n)))+"\n")
                except Exception:
                   print r,n,R,N
                   pass
                
            
def findsiggenepermut(mutfile):
    
    samplelist=[]
    mutdict=defaultdict(list)
    head=0
    #File with all the sample names
    for exp1 in open(mutfile,"rU").xreadlines():
        #print exp1
        lin=exp1.rstrip('\r\n')
        lin=string.split(lin,"\t")
        if len(lin)>2:
            if head==0:
                for i in lin[1:]:
                    i= string.split(i,".")[0]
                    samplelist.append(i)
                head=1
                continue
            else:
                
                for j in range(1,len(lin)):
                    if lin[j]==str(1):
                        mutdict[lin[0]].append(samplelist[j-1])
                        
        else:
            mutdict[lin[1]].append(lin[0])
   # print mutdict           
    return  mutdict

                    
                
if __name__ == '__main__':

    import getopt
  
    mutdict=defaultdict(list)
    
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Insufficient command line flags supplied."
        sys.exit()
    else:
        analysisType = []

        options, remainder = getopt.getopt(sys.argv[1:],'', ['Guidefile=','Reference=','Expand='])
        for opt, arg in options:
            if opt == '--Guidefile': Guidefile=arg
            elif opt == '--Reference':Reference=arg
            elif opt =='--Expand': Expand=arg
            
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
    mutfile=Reference
    #mutfile="/Users/meenakshi/Documents/Mutation_template.txt"
    #mutfile="/Volumes/Pass/Leucegene_Complete/ExpressionInput/round2_sorted/exp.splicing-filteredcor_depleted-filteredNMFsnmf_binary30.txt"
    #Guidefile="/Users/meenakshi/Documents/leucegene/ICGS/Round2_cor_0.6_280default/Clustering-exp.round2_insignificantU2like-Guide1 DDX5&ENSG00000108654&E3.4-E3.9__ENSG0000010-hierarchical_cosine_correlation.txt"          
    header=header_file(mutfile)
    mutdict=findsiggenepermut(mutfile)
    
        #header=header_file(Guidefile)
    Enrichment(Guidefile,mutdict,mutfile,Expand,header)

        
        
        