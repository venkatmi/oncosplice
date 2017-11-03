#!/usr/bin/env python

#!/usr/bin/env python
import traceback
import export
try:
    import math
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",category=UserWarning) ### hides import warnings
        import matplotlib
        matplotlib.use('TkAgg')
        if commandLine==False:
            try: matplotlib.rcParams['backend'] = 'TkAgg'
            except Exception: pass
        try:
            import matplotlib.pyplot as pylab
            import matplotlib.colors as mc
            import matplotlib.mlab as mlab
            from matplotlib import mpl
            from matplotlib.patches import Circle
            from mpl_toolkits.mplot3d import Axes3D
            mpl.rcParams['axes.linewidth'] = 0.5
            mpl.rcParams['pdf.fonttype'] = 42
            mpl.rcParams['font.family'] = 'sans-serif'
            mpl.rcParams['font.sans-serif'] = 'Arial'
        except Exception:
            print 'Matplotlib support not enabled'
        import scipy
        from scipy.linalg import svd
        import scipy.cluster.hierarchy as sch
        import scipy.spatial.distance as dist
        try: import numpy; np = numpy
        except Exception:
            print 'Numpy import error...'
            print traceback.format_exc()
        try:
            import igraph.vendor.texttable
        except ImportError: pass
        try:
            from sklearn.decomposition import PCA, FastICA
        except Exception: pass
        #pylab.ion() # closes Tk window after show - could be nice to include
except Exception:
    print traceback.format_exc()
    pass
import numpy as np
#import pylab as pl
import sys,string
import os
import os.path
from collections import defaultdict
from sklearn.cluster import KMeans
import nimfa
#import statistics

from sklearn.decomposition import NMF

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

def FilterFile(Guidefile,Guidefile_block,PSI):
    if 'Clustering' in Guidefile:
        count=1
        flag=True
        rank_Count=0
        prev=0
    else:
        count=0
    val=[]
    head=0
    for line in open(Guidefile_block,'rU').xreadlines():
        if head >count:
            
            line=line.rstrip('\r\n')
            q= string.split(line,'\t')
            #val.append(q[0])
            if flag:
               
                if int(q[1])==prev:
                    continue
                else:
                    rank_Count+=1
                    prev=int(q[1])
                
        else:
            head+=1
            continue
    head=0
    for line in open(Guidefile,'rU').xreadlines():
        if head >count:
            
            line=line.rstrip('\r\n')
            q= string.split(line,'\t')
            val.append(q[0])
            
                
        else:
            head+=1
            continue
    output_file = PSI[:-4]+'-filtered.txt'
    filterRows(PSI,output_file,filterDB=val)
    
    return output_file,rank_Count


def NMFAnalysis(filename,Rank,turn=0):
    
    X=[]
    header=[]
    head=0
    exportnam=filename[:-4]+'NMFsnmf_versionr'+str(Rank)+'.txt'
    export_res=open(exportnam,"w")
    exportnam1=filename[:-4]+'NMFsnmf_binary'+str(Rank)+'.txt'
    export_res1=open(exportnam1,"w")
    exportnam2=export.findParentDir(filename)+'/round'+str(turn)+'/'+'Metadata'+str(Rank)+'.txt'
    #exportnam2=filename[:-4]+'Metadata'+str(Rank)+'.txt'
    export_res2=export.ExportFile(exportnam2)
    exportnam3=export.findParentDir(filename)+'/round'+str(turn)+'/'+'Annotation'+str(Rank)+'.txt'
    #exportnam3=filename[:-4]+'Annotation'+str(Rank)+'.txt'
    export_res3=export.ExportFile(exportnam3)
    if 'Clustering' in filename:
        count=1
        start=2
    else:
        count=0
        start=1
    print Rank
    for line in open(filename,'rU').xreadlines():
        line=line.rstrip('\r\n')
        q= string.split(line,'\t')
        if head >count:
            val=[]
            val2=[]
            me=0.0
            
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
            #if q[1]==prev:
            X.append(val)
          
        else:
            export_res1.write(line)
            export_res.write(line)
            export_res1.write("\n")
            export_res.write("\n")
            header=q
            head+=1
            continue
    
   
    group=defaultdict(list)
        
    sh=[]
    X=np.array(X)
    print X.shape
    mat=[]
    #mat=X
    mat=zip(*X)
    mat=np.array(mat)
    print mat.shape
    #model = NMF(n_components=15, init='random', random_state=0)
    #W = model.fit_transform(mat)
    nmf = nimfa.Snmf(mat,seed="nndsvd", rank=int(Rank), max_iter=20,n_run=10,track_factor=True)
    nmf_fit = nmf()
    W = nmf_fit.basis()
    W=np.array(W)
    print W.shape
    sh=W.shape
    export_res2.write("uid")
    export_res3.write("uid\tUID\tUID\n")
    for i in range(sh[1]):
        export_res2.write("\t"+'V'+str(i))
    export_res2.write("\n")
    if int(Rank)==2:
        par=1
    else:
        par=2
    #for i in range(sh[1]):
    #    val=W[:,i]
    #    me=np.mean(val)
    #    st=np.std(val)
    #    export_res2.write(header[i+1])
    #    for j in range(sh[0]):
    #        if float(W[i][j])>=float(me+(par*st)):
    #          
    #            export_res2.write("\t"+str(1))
    #        else:
    #            export_res2.write("\t"+str(0))
    #       
    #    export_res2.write("\n")
    W=zip(*W)
    W=np.array(W)
    sh=W.shape
    Z=[]
    for i in range(sh[0]):
        new_val=[]
        val=W[i,:]
        num=sum(i > 0.10 for i in val)
        if num >40 or num <3:
            compstd=True
        else:
            compstd=False
        me=np.mean(val)
        st=np.std(val)
        print 'V'+str(i)
        export_res.write('V'+str(i))
        export_res1.write('V'+str(i))
        export_res3.write('V'+str(i)+"\t"+"Covariate"+"\t"+str(1)+"\n")
        for j in range(sh[1]):
            
            if compstd:   
                if float(W[i][j])>=float(me+(par*st)):
                
                    export_res1.write("\t"+str(1))
                    new_val.append(1)
                else:
                    export_res1.write("\t"+str(0))
                    new_val.append(0)
            else:
                if float(W[i][j])>0.1:
                
                    export_res1.write("\t"+str(1))
                    new_val.append(1)
                else:
                    export_res1.write("\t"+str(0))
                    new_val.append(0)
            export_res.write("\t"+str(W[i][j]))
        Z.append(new_val)
        export_res.write("\n")
        export_res1.write("\n")
    Z=zip(*Z)
    Z=np.array(Z)
    sh=Z.shape
    for i in range(sh[0]):
        val=Z[i,:]
        me=np.mean(val)
        st=np.std(val)
        export_res2.write(header[i+1])
        for j in range(sh[1]):
            export_res2.write("\t"+str(Z[i][j]))
        export_res2.write("\n")
    
    
        
        
   
    return exportnam,exportnam1,exportnam2,exportnam3
                           
                
if __name__ == '__main__':

    import getopt
  
    mutdict=defaultdict(list)
    
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Insufficient command line flags supplied."
        sys.exit()
    else:
        analysisType = []

        options, remainder = getopt.getopt(sys.argv[1:],'', ['Guidefile=','Rank=','PSI='])
        for opt, arg in options:
            if opt == '--Guidefile': Guidefile=arg
            #elif opt == '--Rank':Rank=arg
            elif opt == '--PSI':PSI=arg
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
    #mutfile="/Users/meenakshi/Desktop/Leucegene-data1/Mutation_Annotations.txt"          
    #Guidefile="/Users/meenakshi/Documents/leucegene/ICGS/Round2_cor_0.6_280default/Clustering-exp.round2_insignificantU2like-Guide1 DDX5&ENSG00000108654&E3.4-E3.9__ENSG0000010-hierarchical_cosine_correlation.txt"          
    #inputfile,Rank=FilterFile(Guidefile,PSI)
    inputfile=PSI
    Rank=30
    print Rank
    if Rank>1:
       NMFAnalysis(inputfile,Rank)
    else:
        pass
   