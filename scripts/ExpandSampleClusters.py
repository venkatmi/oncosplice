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

from sklearn import datasets, linear_model
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
import Orderedheatmap
from sklearn import svm
from sklearn.multiclass import OneVsOneClassifier
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import LinearSVC
from sklearn import linear_model
import operator
from collections import OrderedDict
from collections import defaultdict
upd_guides=[]

def FindTopUniqueEvents(Guidefile,psi,Guidedir):
    head=0
    guidekeys=[]
    exportnam=os.path.join(Guidedir,"SplicingeventCount1.txt")
    export_class=open(exportnam,"w")
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
            tempkeys[i].sort(key=operator.itemgetter(1),reverse=False)
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
            upd_guides.append(guidekeys[i][0])
    else:
        omitcluster=1
    export_class.write(psi+"\t"+str(len(unique_clusters[0]))+"\n")
    return omitcluster

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

    for line in open(input_file,'rU').xreadlines():
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
                    for it in range(0,len(counter)):
                        orderlst[counter[it]]=line
                    if logData:
                        line = string.join([values[0]]+map(str,(map(lambda x: math.log(float(x)+1,2),values[1:]))),'\t')+'\n'
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

    for i in filterDB:
        event=string.split(i,"|")[0]
        tempevents.append(event)
    for line in open(input_file,'rU').xreadlines():
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
                filteredevents.append(event)
                for it in range(0,len(counter)):
                    orderlst[counter[it]]=line
                if logData:
                    line = string.join([values[0]]+map(str,(map(lambda x: math.log(float(x)+1,2),values[1:]))),'\t')+'\n'

    for i in range(0,len(tempevents)):
        if i in orderlst:
            export_object.write(orderlst[i])
            if "\n" not in orderlst[i]:
                #print i
                export_object.write("\n") 
    export_object.close()
    tempevents2=[]
    for i in range(len(tempevents)):
        if tempevents[i] in filteredevents:
            tempevents2.append(tempevents[i])

    return tempevents2

def findParentDir(filename):
    filename = string.replace(filename,'//','/')
    filename = string.replace(filename,'\\','/')
    x = string.find(filename[::-1],'/')*-1
    return filename[:x]

def Classify(header,Xobs,output_file,grplst,name,turn):
    count=0
    start=1
    Y=[]
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
            me=np.median(val2)
            for i in range(start,len(q)):
                try:
                    val.append(float(q[i]))
                except Exception:
                    val.append(float(me))
            Y.append(val)
        else:
            head+=1
            continue

    Xobs=zip(*Xobs)
    Xobs=np.array(Xobs)
    Xobs=zip(*Xobs)
    Xobs=np.array(Xobs)
    X=grplst
    X=zip(*X)
    X=np.array(X)
    Y=zip(*Y)
    Y=np.array(Y)

    dire = export.findParentDir(export.findParentDir(export.findParentDir(output_file)[:-1])[:-1])
    output_dir = dire+'SVMOutputs'
    if os.path.exists(output_dir)==False:
        export.createExportFolder(output_dir)

    exportnam1=output_dir+'/round'+str(turn)+'SVC_decision_func.txt'
    export_class1=open(exportnam1,"w")
    exportnam2=output_dir+'/round'+str(turn)+'SVC_Results.txt'
    export_class2=open(exportnam2,"w")
    regr = LinearSVC()
    regr.fit(Xobs,X[:,0])
    q=regr.predict(Y)
    count=1

    if len(X[:,0])>2:
        prob_=regr.fit(Xobs,X[:,0]).decision_function(Y)
        export_class1.write("uid")
        export_class2.write("uid")
        for ni in name:
            sub=string.split(ni,"_")[0]
            export_class1.write("\t"+"R"+str(turn)+"-"+sub)
            export_class2.write("\t"+"R"+str(turn)+"-"+sub)
        export_class1.write("\n")
        export_class2.write("\n")

        for iq in range(0,len(header)-1):
            export_class1.write(header[iq+1])
            export_class2.write(header[iq+1])
            for jq in range(0,len(X[:,0])):
                export_class1.write("\t"+str(prob_[iq][jq]))
                if prob_[iq][jq]>0:
                    export_class2.write("\t"+str(1))
                else:
                    export_class2.write("\t"+str(0))
            export_class1.write("\n")
            export_class2.write("\n")
    else:
        prob_=regr.fit(Xobs,X[:,0]).decision_function(Y)
        export_class1.write("uid"+"\t")
        export_class2.write("uid"+"\t")
        export_class1.write("group")
        export_class2.write("R"+str(turn)+"-V1"+"\t"+"R"+str(turn)+"-V2")
        export_class1.write("\n")
        export_class2.write("\n")

        for iq in range(0,len(header)-1):
            export_class1.write(header[iq+1])
            export_class2.write(header[iq+1])
            export_class1.write("\t"+str(prob_[iq]))
            if prob_[iq]>0.5:
                export_class2.write("\t"+str(1)+"\t"+str(0))
            else:
                if prob_[iq]<-0.5:  
                    export_class2.write("\t"+str(0)+"\t"+str(1))
                else:
                    export_class2.write("\t"+str(0)+"\t"+str(0))
            export_class1.write("\n")
            export_class2.write("\n")
    export_class2.close() 
    Orderedheatmap.Classify(exportnam2)
    
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

def avg(array):
    total = sum(map(float, array))
    average = total/len(array)
    return average

def TrainDataGeneration(output_file,NMF_annot,name):
    head=0
    groups=[1,2]
    matrix=defaultdict(list)
    compared_groups={}
    train=[]
    
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
                    try:mapping[0].append(header[i])
                    except Exception: mapping[0]=[header[i]]
            head2=0

            eventname=[]
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
                    gvalues_list=[]
                    for i in group_db[1]:
                        try:
                            x=float(lin2[i])
                            gvalues_list.append(x)
                        except Exception:
                            pass
                    try:  
                        matrix[lin[0]].append(avg(gvalues_list))
                        eventname.append(key)
                    except Exception:
                        matrix[lin[0]].append(float(0))
                        eventname.append(key)
                    
    for j in range(0,len(name)):
        for key in matrix:
            key1=key+"_vs"
            key2="vs_"+key+".txt"
            
            if key1 in name[j] or key2 in name[j]:
                train.append(matrix[key])
    return train

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
    counter=1

    for filename in os.listdir(Guidedir):
        if filename.startswith("PSI."):
            Guidefile=os.path.join(Guidedir, filename)
            psi=string.replace(filename,"PSI.","")
            PSIfile=os.path.join(PSIdir, psi)
            print Guidefile,PSIfile
  
            omitcluster=FindTopUniqueEvents(Guidefile,psi,Guidedir)
            print omitcluster
            if omitcluster==0:
                group.append(counter)
                name.append(psi)
                counter+=1
        
    output_file=PSI[:-4]+"-filtered.txt"  

    print len(upd_guides)
    filterRows(PSI,output_file,filterDB=upd_guides,logData=False)
    header=header_file(output_file)
    
    train=TrainDataGeneration(output_file,NMF_annot,name)
    grplst.append(group)
    print grplst
    Classify(header,train,output_file,grplst,name)