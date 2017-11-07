#!/usr/bin/env python

import numpy as np
import pylab as pl
import sys,string
import os
import os.path
import scipy
import operator
from collections import OrderedDict
from collections import defaultdict
from operator import itemgetter

def Classify(filename):
    count=0
    start=1
    orderdict=OrderedDict()
    countdict=OrderedDict()
  
    countlst=[]
    
    Y=[]
    head=0
    rownames=[]
    colnames=[]
    q=[]
    Z=[]
    output_file=filename[:-4]+"-filtered.txt"
    export_object = open(output_file,'w')
    for line in open(filename,'rU').xreadlines():
        if head >0:
            val=[]
            counter2=0
            val2=[]
            me=0.0
            
            line=line.rstrip('\r\n')
            
            q= string.split(line,'\t')
           # rownames.append(q[0])
            orderdict[q[0]]=[q[0],]
            for i in range(start,len(q)):
                try:
                    val2.append(float(q[i]))
                    try:
                        orderdict[q[0]].append(float(q[i]))
                    except Exception:
                        orderdict[q[0]]=[float(q[i]),]
                    try:
                        countdict[i].append(float(q[i]))
                    except Exception:
                        countdict[i]=[float(q[i]),]
                except Exception:
                    continue
            
            
            count+=1
        else:
            #export_object.write(line)
            head=1
            line=line.rstrip('\r\n')
            
            q= string.split(line,'\t')
            header=q
            continue
    
    for i in countdict:
       
        countlst.append(sum(countdict[i]))
    print countlst
    
    B=sorted(range(len(countlst)),key=lambda x:countlst[x],reverse=True)
    C=sorted(range(len(countlst)),key=lambda x:B[x])
    print C
    qu=0
    for i in orderdict.keys():
        Y.append(orderdict[i])
        qu+=1
        print qu
    
    for i in range(0,len(C)):
        jk= C.index(i)+1
        print jk
        print Y[jk]
        Y=sorted(Y,key=itemgetter(jk))
        
        #orderdict=OrderedDict(sorted(orderdict,key=itemgetter(jk)))
        #colnames.append(header[C.index(i)+1])
    
    Y=np.array(Y)
    Y=zip(*Y)
    Y=np.array(Y)
    Z.append(Y[0,:])
    for i in range(0,len(C)):
        jk= C.index(i)+1
        Z.append(Y[jk,:])
    Z=np.array(Z)
    q= Z.shape
    print Z.shape
    print len(header)
    #print header[24]
    #print q[0]
    #print q[1]
    #for i in range(0,len(C)):
    #    Y=sorted(Y,key=itemgetter(i))
    export_object.write("uid")
    #for i in colnames:
    #    export_object.write("\t"+i)
    #export_object.write("\n")
        
    for i in range(q[1]):
        export_object.write("\t"+Z[0][i])
    export_object.write("\n")
    for ij in range(1,q[0]):
        jk= C.index(ij-1)+1
        export_object.write(header[jk])
        for jq in range(0,q[1]):
            export_object.write("\t"+str(Z[ij][jq]))
        export_object.write("\n")
        
    
    
    
    
            
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

        options, remainder = getopt.getopt(sys.argv[1:],'', ['Guidedir='])
        for opt, arg in options:
            if opt == '--Guidedir': Guidedir=arg
           
           
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
                
    Classify(Guidedir)