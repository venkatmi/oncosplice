#!/usr/bin/env python

import sys, string, os

import ExpandClusters_SVM as ExpandSampleClusters


import numpy as np
import export
upd_guides=[]
import operator
from collections import OrderedDict
from collections import defaultdict


def FindTopUniqueEvents(Guidefile,psi,Guidedir,adjp):
    head=0
    guidekeys=[]
    
    #commonkeys=[]
    tempkeys={}
    global upd_guides
    global train
    omitcluster=0
    
    unique_clusters={}

    for line in open(Guidefile,'rU').xreadlines():
      
        if head==0:
            head=1
            continue
        else:
            line1=line.rstrip('\r\n')
            q= string.split(line1,'\t')
            
            if adjp:
                if (abs(float(q[9]))>0.1 and float(q[11])<0.05) :
                    try:
                        tempkeys[q[4]].append([q[0],float(q[11]),q[12]])
                    except KeyError:
                        tempkeys[q[4]]=[[q[0],float(q[11]),q[12]],]
            else:
                if (abs(float(q[9]))>0.1 and float(q[10])<0.05) :
                    try:
                        tempkeys[q[4]].append([q[0],float(q[11]),q[12]])
                    except KeyError:
                        tempkeys[q[4]]=[[q[0],float(q[11]),q[12]],]
    for i in tempkeys:
        
        if len(tempkeys[i])>1:
            #print tempkeys[i]
            tempkeys[i].sort(key=operator.itemgetter(1),reverse=False)
            #print tempkeys[i][0]
            try:
                unique_clusters[0].append(tempkeys[i][0])
            except KeyError:
                unique_clusters[0]=[tempkeys[i][0],]
          
        else:
            try:
                unique_clusters[0].append(tempkeys[i][0])
            except KeyError:
                unique_clusters[0]=[tempkeys[i][0],]
    
    try:
        
        if len(unique_clusters[0])>100:     
            unique_clusters[0].sort(key=operator.itemgetter(1))
  
            guidekeys=unique_clusters[0]
            for i in range(0,len(guidekeys)):
       
                    upd_guides.append(guidekeys[i][0])
           
            
        else:
            omitcluster=1
    except Exception:
        omitcluster=1
    return omitcluster

def normalization(filename):
    output_file=filename[:-4]+"-normalized.txt"
    exportnam=open(output_file,"w")
    head=0
    for line in open(filename,'rU').xreadlines():
        if head==0:
            exportnam.write(line)
            head=1
            continue
        else:
            val2=[]
            me=0.0
            line=line.rstrip('\r\n')
            q= string.split(line,'\t')
            for i in range(1,len(q)):
                try:
                    val2.append(float(q[i]))
                except Exception:
                    continue
            me=np.median(val2)
            exportnam.write(q[0])
            for i in range(1,len(q)):
                try:val=float(q[i])-me
                except Exception: val= me
                exportnam.write("\t"+str(val))
            exportnam.write("\n")
            #if q[1]==prev:
    exportnam.close()       
    return output_file


def SupervisedAnalyses(Training,Test,Multigroup,group,diffevents,adjp,centroid,normalize):
    
        species="Hs"
        counter=1
        Guidedir=diffevents 
        BinarizedOutput=group
        global upd_guides
        upd_guides=[]
        name=[]
        
        grplst=[]
        filteredevents=[]
       
        for filename in os.listdir(Guidedir):
            
            counter=1
            if filename.startswith("PSI."):
                name=[]
                upd_guides=[]
                Guidefile=os.path.join(Guidedir, filename)
                psi=string.replace(filename,"PSI.","")
               
                omitcluster=FindTopUniqueEvents(Guidefile,psi,Guidedir,adjp) 
                if omitcluster==0:
                    
                    name.append(psi)
                    counter+=1
                InputFile2=Training
                InputFile=Test
                output_file=InputFile[:-4]+"-filtered.txt"
                output_file1=InputFile2[:-4]+"-filtered.txt" 
            
                filteredevents=ExpandSampleClusters.filterRows_data(InputFile,output_file,filterDB=upd_guides,logData=False)
                ExpandSampleClusters.filterRows_data(InputFile2,output_file1,filterDB=filteredevents,logData=False)
                
                if normalize==True:
                    output_file=normalization(output_file)
                    output_file1=normalization(output_file1)
                header=ExpandSampleClusters.header_file(output_file)
                header1=ExpandSampleClusters.header_file(output_file1)
         
                train,train1=ExpandSampleClusters.TrainDataGeneration(output_file1,BinarizedOutput,header,name,Multigroup)
                grplst=[1,0]
                ExpandSampleClusters.Classify(header,train,train1,output_file,grplst,name)

if __name__ == '__main__':
    import getopt
    adjp=False
    centroid=True
    Multigroup=False
    normalize=False
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Insufficient command line flags supplied."
        sys.exit()
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['Training=','Test=','Multigroup=','group=','diffevents=','adjp=','centroid=','normalize='])
        for opt, arg in options:
            if opt == '--Training': Training=arg
            elif opt=='--Test': Test=arg
            elif opt=='--Multigroup': Multigroup=arg
            elif opt=='--group': group=arg
            elif opt=='--diffevents': diffevents=arg
            elif opt=='--adjp': adjp=arg
            elif opt=='--centroid':centroid=arg
            elif opt=='--normalize':normalize=arg
            
    SupervisedAnalyses(Training,Test,Multigroup,group,diffevents,adjp,centroid,normalize)
  
    print "completed"