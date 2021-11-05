#!/usr/bin/env python

import numpy as np
import pylab as pl
import sys,string
import os
import os.path
from collections import defaultdict
from sklearn.cluster import KMeans

try:from stats_scripts import statistics
except Exception: import statistics
import random
import UI
import export; reload(export)
import re
try:from stats_scripts import fishers_exact_test
except Exception:import fishers_exact_test
import traceback
import warnings
import math
import export

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
            z = (r - (n*(R/N)))/math.sqrt(n*(R/N)*(1-(R/N))*(1-((n-1)/(N-1))))
            print r, n, R, N, z;sys.exit()
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
        print ft.two_tail_p()
        return ft.two_tail_p(),z
    

            
def returnSamplesInMetaData(fname, delimiter=None,metaDataMatrixFormat=False):
    """ Returns the samples present based on the two possible metadata file input formats
    (three column groups file or binary matrix of samples by mutations) """
    head=0
    header=[]
    newheader=[]
    with open(fname, 'rU') as fin:
        for line in fin:
            #print line
            line = line.rstrip(os.linesep)
            line = string.replace(line,'"','')
            header=string.split(line,'\t')
        
            if metaDataMatrixFormat==True:
                if head==0:
                    
                    for i in range(1,len(header)):
                        iq=header[i]
                        #iq=string.split(header[i],".")[0]
                        newheader.append(iq)
                    head=1
                else:break
            else:
                if len(header)<3 or metaDataMatrixFormat==False:
                    if header[0] not in newheader:
                        newheader.append(header[0]) 
    #print len(newheader)           
    return newheader

def Enrichment(Inputfile,mutdict,mutfile,metaDataMatrixFormat,header):
    import collections
    import mappfinder
    X=defaultdict(list)
    prev=""
    head=0
    group=defaultdict(list)
    enrichdict=defaultdict(float)
    mut=export.findFilename(mutfile)
    dire=export.findParentDir(Inputfile)
    output_dir = dire+'MutationEnrichment'
    print output_dir
    export.createExportFolder(output_dir)
    number_of_samples = 0
    
    ### All enrichment results
    exportnam=output_dir+'/Enrichment_Results.txt'
    export_enrich=open(exportnam,"w")
    
    ### Selected Enrichment results based on p-value, sensitivity and specificity for association with cluster names
    exportnam=output_dir+'/Enrichment_tophits.txt'
    export_hit=open(exportnam,"w")
   
    header = "Mutations"+"\t"+"Cluster"+"\t"+"r"+"\t"+"R"+"\t"+"n"+"\t"+"Sensitivity"+"\t"+"Specificity"+"\t"+"z-score"+"\t"+"Fisher exact test"+"\t"+"adjp value"+"\n"
    export_enrich.write(header)
    export_hit.write(header)
    header2=returnSamplesInMetaData(Inputfile,metaDataMatrixFormat=True)
    print header2
    for line in open(Inputfile,'rU').xreadlines():
        if head > 0:
            number_of_samples+=1
            line=line.rstrip('\r\n')
            q = string.split(line,'\t')
            for i in range(1,len(q)):
                if q[i]==str(1):
                    #group[q[0]].append(header2[i-1])
                    group[header2[i-1]].append(q[0]) ### [Cluster] = [full_sample_ID]
        else:
            head+=1
            continue
   
    print 'Number of patient samples in dataset =',number_of_samples
    total_Scores={}
    for kiy in mutdict:
        if kiy =="MDP":
            print mutdict[kiy]
        groupdict={}
        remaining=[]
        remaining=list(set(header) - set(mutdict[kiy]))
        groupdict[1]=mutdict[kiy]
        groupdict[2]=remaining
        #export_enrich1.write(kiy)
        for key2 in group:
            r=float(len(list(set(group[key2])))-len(list(set(group[key2]) - set(mutdict[kiy]))))
            n=float(len(group[key2]))
            R=float(len(set(mutdict[kiy])))
            N=float(number_of_samples)
            if r==0 or key2=="1" or R==1.0:
                #print kiy,key2,r,n,R,N
                pval=float(1)
                z=float(0)
                null_z = 0.000
                zsd = mappfinder.ZScoreData(key2,r,R,z,null_z,n)
                zsd.SetP(pval)
            else:
                try: z = Zscore(r,n,N,R)
                except: z=0
                ### Calculate a Z-score assuming zero matching entries
                try: null_z = Zscore(0,n,N,R)
                except Exception: null_z = 0.000
               
                try:
                    pval = mappfinder.FishersExactTest(r,n,R,N)
                    zsd = mappfinder.ZScoreData(key2,r,R,z,null_z,n)
                    zsd.SetP(pval)
                except Exception:
                    pval=1.0
                    zsd = mappfinder.ZScoreData(key2,r,R,z,null_z,n)
                    zsd.SetP(pval)
                    #pass
                
            if kiy in total_Scores:
                    signature_db = total_Scores[kiy]
                    signature_db[key2]=zsd ### Necessary format for the permutation function
            else:
                    signature_db={key2:zsd}
                    total_Scores[kiy] = signature_db
    sorted_results=[]
    mutlabels={}
    for kiy in total_Scores:
        signature_db = total_Scores[kiy]
        ### Updates the adjusted p-value instances
        mappfinder.adjustPermuteStats(signature_db)
        for signature in signature_db:
            zsd = signature_db[signature]
            results = [kiy,signature,zsd.Changed(),zsd.Measured(),zsd.InPathway(),str(float(zsd.PercentChanged())/100.0),str(float(float(zsd.Changed())/float(zsd.InPathway()))), zsd.ZScore(), zsd.PermuteP(), zsd.AdjP()] #string.join(zsd.AssociatedIDs(),'|')
            sorted_results.append([signature,-1*float(zsd.ZScore()),results])
    sorted_results.sort() ### Sort z-score
    
    prev=""
    for (sig,p,values) in sorted_results:
        if sig!=prev:
            flag=True
            export_hit.write(string.join(values,'\t')+'\n')
        if flag:
            ### Update the cluster label to include the top enriched term meeting, sensitivity and specificity cutoffs
            #print values[5],values[6],values[6],values[2]; sys.exit()
            if (float(values[5])>=0.2 and float(values[6])>=0.2 and float(values[7])>=1.95 and float(values[2])>=2):
                clusterID = values[1]
                topEnrichedTerm=values[0]
                mutlabels[clusterID]=clusterID+' ('+topEnrichedTerm+')'
                flag=False
                export_hit.write(string.join(values,'\t')+'\n')
        export_enrich.write(string.join(values,'\t')+'\n')
        prev=sig
    if len(sorted_results)==0:
            export_enrich.write(string.join([splicing_factor,'NONE','NONE','NONE','NONE','NONE','NONE'],'\t')+'\n')
    export_enrich.close()

    return mutlabels
            
def findsiggenepermut(mutfile,valid_filenames=None):
    samplelist=[]
    mutdict=defaultdict(list)
    head=0
    #File with all the sample names
    for exp1 in open(mutfile,"rU").xreadlines():
        #print exp1
        lin=exp1.rstrip('\r\n')
        lin = string.replace(lin,'"','')
        lin=string.split(lin,"\t")
        if len(lin)>3:
            if head==0:
                for i in lin[1:]:
                    samplelist.append(i)
                head=1
                continue
            else:
                for j in range(1,len(lin)):
                    if lin[j]==str(1):
                        sampleName = samplelist[j-1]
                        if valid_filenames != None:
                            for fullName in valid_filenames:
                                if sampleName in fullName:
                                    sampleName = fullName
                        mutdict[lin[0]].append(fullName)
        else:
            sampleName = lin[0]
            if valid_filenames != None:
                for fullName in valid_filenames:
                    if sampleName in fullName:
                        sampleName = fullName
            mutdict[lin[-1]].append(sampleName)
    return  mutdict

def Zscore(r,n,N,R):
    """where N is the total number of events measured: 
    R is the total number of events meeting the criterion:
    n is the total number of events in this specific reference gene-set: 
    r is the number of events meeting the criterion in the examined reference gene-set: """
    N=float(N) ### This bring all other values into float space

    z = (r - n*(R/N))/math.sqrt(n*(R/N)*(1-(R/N))*(1-((n-1)/(N-1))))
    return z               
                
if __name__ == '__main__':
    """ Enrichment analysis for user-supplied metadata, including mutations and clincial characteristics """
    metaDataMatrixFormat = False
    import getopt
    mutdict=defaultdict(list)
    
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Insufficient command line flags supplied."
        sys.exit()
    else:
        analysisType = []

        options, remainder = getopt.getopt(sys.argv[1:],'', ['Inputfile=','Reference=','metaDataMatrixFormat=', 'i=', 'input=', 'r=', 'reference=', 'expand='])
        for opt, arg in options:
            if opt == '--Inputfile' or opt == '--i' or opt == '--input': Inputfile=arg
            elif opt == '--Reference' or opt == '--r' or opt == '--reference': Reference=arg
            elif opt =='--metaDataMatrixFormat' or opt =='--expand':
                if string.lower(arg) == 'true' or string.lower(arg) == 'yes':
                    metaDataMatrixFormat = True
                else:
                    metaDataMatrixFormat = False
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
    mutfile=Reference
    
    header=returnSamplesInMetaData(mutfile,metaDataMatrixFormat=metaDataMatrixFormat)
    mutdict=findsiggenepermut(mutfile)
    Enrichment(Inputfile,mutdict,mutfile,metaDataMatrixFormat,header)

        
        
        