#!/usr/bin/env python

#!/usr/bin/env python

import numpy as np
import pylab as pl
import sys,string
import os
import os.path
import scipy
import sampleIndexSelection
import matplotlib.pyplot as plt
import export
import os.path
from os import path
from sklearn import datasets, linear_model
from sklearn.preprocessing import StandardScaler

from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC

from sklearn import svm
from sklearn.multiclass import OneVsOneClassifier
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import LinearSVC
from sklearn import linear_model
import operator
from collections import OrderedDict
from collections import defaultdict

upd_guides=[]

#upd_guides.append("uid")
def FindTopUniqueEvents(Guidefile,psi,Guidedir):
    head=0
    guidekeys=[]
    exportnam=os.path.join(Guidedir,"SplicingeventCount1.txt")
    export_class=open(exportnam,"w")
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
            
            if abs(float(q[8]))>0.15:
                try:
                    tempkeys[q[2]].append([q[0],float(q[10]),q[11]])
                except KeyError:
                    tempkeys[q[2]]=[[q[0],float(q[10]),q[11]],]
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
          
    unique_clusters[0].sort(key=operator.itemgetter(1))
  
    if len(unique_clusters[0])>100:
        guidekeys=unique_clusters[0]
        for i in range(0,len(guidekeys)):
            
            #upd_guides[i]=[upd_guides[i][3],upd_guides[i][4]]
            upd_guides.append(guidekeys[i][0])
    else:
        omitcluster=1
    export_class.write(psi+"\t"+str(len(unique_clusters[0]))+"\n")
    print len(upd_guides)
    return omitcluster
    #return upd_guides,train

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def filterRows(input_file,output_file,filterDB=None,logData=False):
    orderlst={}
    counter=[]
    export_object = open(output_file,'w')
    firstLine = True
    Flag=0;
    print len(filterDB)
    #for i in filterDB:
    for line in open(input_file,'rU').xreadlines():
        #for i in filterDB:
            flag1=0
            
            data = cleanUpLine(line)
           
            values = string.split(data,'\t')
          
            
            if firstLine:
                firstLine = False
                if Flag==0:
                    export_object.write(line)
            else:
                if values[0] in filterDB:
                    counter=[index for index, value in enumerate(filterDB) if value == values[0]]
                        #print counter
                    for it in range(0,len(counter)):
                        orderlst[counter[it]]=line
                    if logData:
                        line = string.join([values[0]]+map(str,(map(lambda x: math.log(float(x)+1,2),values[1:]))),'\t')+'\n'
                        #export_object.write(line)
                        #firstLine=True
                       # Flag=1;
               
                    
                #else:
                   # max_val = max(map(float,values[1:]))
                #min_val = min(map(float,values[1:]))
                #if max_val>0.1:

                     #   export_object.write(line)
    try:
        for i in range(0,len(orderlst)):
            export_object.write(orderlst[i])
    except Exception:
        print i,filterDB[i]
    
    
         
    export_object.close()
    print 'Filtered rows printed to:',output_file
    
def filterRows_data(input_file,output_file,filterDB=None,logData=False):
    filteredevents=[]
    tempevents=[]
    orderlst={}
    counter=[]
    export_object = open(output_file,'w')
    firstLine = True
    Flag=0;
    print len(filterDB)
    for i in filterDB:
        event=string.split(i,"|")[0]
        tempevents.append(event)
    for line in open(input_file,'rU').xreadlines():
        #for i in filterDB:
            flag1=0
            
            data = cleanUpLine(line)
           
            values = string.split(data,'\t')
            event=string.split(values[0],"|")[0]
                  
            
            if firstLine:
                firstLine = False
                if Flag==0:
                    export_object.write(line)
            else:
                if event in tempevents:
                    counter=[index for index, value in enumerate(tempevents) if value == event]
                    #print counter
                    filteredevents.append(event)
                    for it in range(0,len(counter)):
                        orderlst[counter[it]]=line
                    if logData:
                        line = string.join([values[0]]+map(str,(map(lambda x: math.log(float(x)+1,2),values[1:]))),'\t')+'\n'
                        #export_object.write(line)
                        #firstLine=True
                       # Flag=1;
               
                    
                #else:
                   # max_val = max(map(float,values[1:]))
                #min_val = min(map(float,values[1:]))
                #if max_val>0.1:
    
                     #   export_object.write(line)
   
    for i in range(0,len(tempevents)):
        if i in orderlst:
            export_object.write(orderlst[i])

    
    
         
    export_object.close()
    tempevents2=[]
    print 'Filtered rows printed to:',output_file
    for i in range(len(tempevents)):
        if tempevents[i] in filteredevents:
            tempevents2.append(tempevents[i])
    
        
    return tempevents2

        
def header_file(fname, delimiter=None):
    head=0
    header=[]
    new_head=[]
    with open(fname, 'rU') as fin:
        for line in fin:
            if head==0:
                line = line.rstrip(os.linesep)
                header=string.split(line,'\t')
                for i in header:
                    if ":" in i:
                        
                        i=string.split(i,":")[1]
                    new_head.append(i)
                        
                head=1
            else:break
   
    return new_head

def Classify(header,Xobs1,train2,output_file,grplst,name):
    count=0
    start=1
    Y=[]
    header1=[]
    head=0
    for line in open(output_file,'rU').xreadlines():
        if head >count:
            val=[]
            counter2=0
            val2=[]
            me=0.0
            line=line.rstrip('\r\n')
            q= string.split(line,'\t')
            for i in range(start,len(q)):
                try:
                    val2.append(float(q[i]))
                except Exception:
                    continue
            me=np.mean(val2)
            for i in range(start,len(q)):
                try:
                    val.append(float(q[i]))
                except Exception:
                    val.append(float(me))
            #if q[1]==prev:
            Y.append(val)
        
        else:
            header=line
            line=line.rstrip('\r\n')
            header1=string.split(line,"\t")
            head+=1
            continue
    #print Xobs
    #exportnam=output_file[:-4]+'KNN_one.txt'
    #export_class_knn=open(exportnam,"a")
    #export_class_knn.write(header)
    exportnam=output_file[:-4]+'_SVMClasses.txt'
    if path.exists(exportnam):
        export_class_svc=open(exportnam,"a")
        exportnam1=output_file[:-4]+'_SVMDecisionFunc.txt'
        export_class1=open(exportnam1,"a")
    else:
        export_class_svc=open(exportnam,"a")
        exportnam1=output_file[:-4]+'_SVMDecisionFunc.txt'
        export_class1=open(exportnam1,"a")
        export_class_svc.write(header)
    
        for iqt in range(0,len(header1)-1):
            export_class1.write("\t"+header1[iqt+1])
        export_class1.write("\n")
 
    Y=zip(*Y)
    Y=np.array(Y)

  
    for iq in range(0,len(name)):
        Xobs=[Xobs1[iq],train2[iq]]
    
        Xobs=zip(*Xobs)
   
        Xobs=np.array(Xobs)
        Xobs=zip(*Xobs)

        Xobs=np.array(Xobs)
        X=[[1,2]]
        X=zip(*X)
        X=np.array(X)
       
        
        regr = LinearSVC()
        regr.fit(Xobs,X[:,0])
        q=regr.predict(Y)

        count=1
        export_class_svc.write(name[iq])
        for i in q:
        
            export_class_svc.write("\t"+str(i))
            count+=1
        export_class_svc.write("\n")
        prob_=regr.fit(Xobs,X[:,0]).decision_function(Y)
      
        export_class1.write("\n")
     
     
       
        export_class1.write(name[iq])
        for iqt in range(0,len(header1)-1):
           # export_class1.write(header1[iqt+1])
           
            #for jq in range(0,len(X[:,0])):
            export_class1.write("\t"+str(prob_[iqt]))
            
        export_class1.write("\n")
        
       
       
        
    

def avg(array):
    total = sum(map(float, array))
    average = total/len(array)
    return average

def TrainDataGeneration(output_file,NMF_annot,header,name,Multigroup):
    head=0
    groups=[1,2]
    matrix=defaultdict(list)
    matrix2=defaultdict(list)
    compared_groups={}
    train=[]
    train2=[]
    keye="None"
    keye2="None"
    
    
    if Multigroup==True:
     for exp1 in open(NMF_annot,"rU").xreadlines():
        lin=exp1.rstrip('\r\n')
        lin=string.split(lin,"\t")
        mapping={}
        if head==0:
            header=lin
            head=1
            continue
        else:
            for i in range(1,len(lin)):
                if lin[i]=='1':
                    try:mapping[1].append(header[i])
                    except Exception: mapping[1]=[header[i]]
                else:
                    #if lin[i]=='0':
                        #try:mapping[3].append(header[i])
                        #except Exception: mapping[3]=[header[i]]
                        try:mapping[0].append(header[i])
                        except Exception: mapping[0]=[header[i]]
            head2=0
            
           
            for exp2 in open(output_file,"rU").xreadlines():
                    lin2=exp2.rstrip('\r\n')
                    lin2=string.split(lin2,"\t")
                    if head2==0:
                        group_db={}
                        index=0
                        try:
                            if len(mapping[1])>0 and len(mapping[0])>0:
                                for i in lin2[1:]:
                                    if i in mapping[1]:
                                        try: group_db[1].append(index)
                                        except Exception: group_db[1] = [index]
                                    else:
                                        if i in mapping[0]:
                                            try: group_db[2].append(index)
                                            except Exception: group_db[2] = [index]
                                    index+=1
                        except Exception:
                            break
                                
                            
                        head2=1
                        continue
                    else:
                        key = lin2[0]
                        lin2=lin2[1:]
                        grouped_floats=[]
                        associated_groups=[]
            
                        ### string values
                        gvalues_list=[]
                        gvalues_list2=[]
                        for i in group_db[1]:
                                try:
                                    
                                    x=float(lin2[i])
                                    gvalues_list.append(x)
                                            
                                          #  math.log((float(values[i])+1,2))
                                
                                    
                                except Exception:
                                    #try: gvalues_list.append('') ### Thus are missing values
                                    #except Exception: pass
                                    pass
                        try:
                        
                            matrix[lin[0]].append(avg(gvalues_list))
                        except Exception:
                            matrix[lin[0]].append(float(0))
                        gval_list=[]
                        for i in group_db[2]:
                                try:
                                    
                                    x=float(lin2[i])
                                    gval_list.append(x)
                                            
                                          #  math.log((float(values[i])+1,2))
                                
                                    
                                except Exception:
                                    #try: gvalues_list.append('') ### Thus are missing values
                                    #except Exception: pass
                                    pass
                        try:
                        
                            matrix2[lin[0]].append(avg(gval_list))
                        except Exception:
                            matrix2[lin[0]].append(float(0))
                            
    else:
     mapping={} 
     for exp1 in open(NMF_annot,"rU").xreadlines():
        lin=exp1.rstrip('\r\n')
        lin=string.split(lin,"\t")
        if lin[1]=='1':
            try:mapping[1].append(lin[0])
            except Exception: mapping[1]=[lin[0]]
            ky=lin[2]
        else:
            if lin[1]=='2':
                try:mapping[0].append(lin[0])
                except Exception: mapping[0]=[lin[0]]
            ky2=lin[2]
     keye=ky+"_vs_"+ky2
     keye2=ky2+"_vs_"+ky      
            
     head2=0
     for exp2 in open(output_file,"rU").xreadlines():
        lin2=exp2.rstrip('\r\n')
        lin2=string.split(lin2,"\t")
        if head2==0:
            group_db={}
            index=0
            try:
                if len(mapping[1])>0 and len(mapping[0])>0:
                    for i in lin2[1:]:
                        if i in mapping[1]:
                            try: group_db[1].append(index)
                            except Exception: group_db[1] = [index]
                        else:
                            if i in mapping[0]:
                                try: group_db[2].append(index)
                                except Exception: group_db[2] = [index]
                        index+=1
            except Exception:
                break
                    
                
            head2=1
            continue
        else:
            key = lin2[0]
            lin2=lin2[1:]
            grouped_floats=[]
            associated_groups=[]

            ### string values
            gvalues_list=[]
            gvalues_list2=[]
            for i in group_db[1]:
                    try:
                        
                        x=float(lin2[i])
                        gvalues_list.append(x)
                                
                              #  math.log((float(values[i])+1,2))
                    
                        
                    except Exception:
                        #try: gvalues_list.append('') ### Thus are missing values
                        #except Exception: pass
                        pass
            try:
            
                matrix[keye].append(avg(gvalues_list))
            except Exception:
                matrix[keye].append(float(0))
            gval_list=[]
            for i in group_db[2]:
                    try:
                        
                        x=float(lin2[i])
                        gval_list.append(x)
                                
                              #  math.log((float(values[i])+1,2))
                    
                        
                    except Exception:
                        #try: gvalues_list.append('') ### Thus are missing values
                        #except Exception: pass
                        pass
            try:
            
                matrix2[keye].append(avg(gval_list))
            except Exception:
                matrix2[keye].append(float(0))
                       
                        #
   
    for j in range(0,len(name)):
       
        for key in matrix:
           
            key1=key+"_vs"
            if key1 in name[j] or keye in name[j] or keye2 in name[j]:
               
                
                train.append(matrix[key])
                
        for key in matrix2:
            key1=key+"_vs"
            if key1 in name[j] or keye in name[j] or keye2 in name[j]:
                
                train2.append(matrix2[key])
              
    #print train,train2
    return train,train2



if __name__ == '__main__':

    import getopt
    group=[]
    grplst=[]
    name=[]
    matrix={}
    compared_groups={}

    
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Insufficient command line flags supplied."
        sys.exit()
    else:
        analysisType = []

        options, remainder = getopt.getopt(sys.argv[1:],'', ['Guidedir=','PSIdir=','PSI=','NMF_annot='])
        for opt, arg in options:
            if opt == '--Guidedir': Guidedir=arg
            elif opt =='--PSIdir':PSIdir=arg
            elif opt =='--PSI':PSI=arg
            elif opt =='--NMF_annot':NMF_annot=arg
           
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
    #commonkeys=[]
    counter=1
#filename="/Users/meenakshi/Documents/leucegene/ICGS/Clustering-exp.Hs_RNASeq_top_alt_junctions367-Leucegene-75p_no149-Guide1 TRAK1&ENSG00000182606&I1.1_42075542-E2.1__E-hierarchical_cosine_correlation.txt"          
#PSIfile="/Users/meenakshi/Documents/leucegene/ExpressionInput/exp.Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation-367-Leucegene-75p-unique-filtered-filtered.txt"
#keylabel="/Users/meenakshi/Documents/leucegene/ExpressionInput/exp.round2_glmfilteredKmeans_label.txt"
    for filename in os.listdir(Guidedir):
        if filename.startswith("PSI."):
            Guidefile=os.path.join(Guidedir, filename)
            psi=string.replace(filename,"PSI.","")
            PSIfile=os.path.join(PSIdir, psi)
            print Guidefile,PSIfile
  
            #output_file=PSIfile[:-4]+"-filtered.txt"
    #sampleIndexSelection.filterFile(PSIfile,output_file,header)
            omitcluster=FindTopUniqueEvents(Guidefile,psi,Guidedir)
            print omitcluster
            if omitcluster==0:
                group.append(counter)
                name.append(psi)
                counter+=1
        
    output_file=PSI[:-4]+"-filtered.txt"  
            #print guidekey
    print len(upd_guides)
    filterRows(PSI,output_file,filterDB=upd_guides,logData=False)
    header=header_file(output_file)
    
    train=TrainDataGeneration(output_file,NMF_annot,name)
    grplst.append(group)
    print grplst
    Classify(header,train,output_file,grplst,name)