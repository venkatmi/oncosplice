#!/usr/bin/env python

#Copyright 2017 Cincinnati Children's Hospital Medical Center, Research Foundation
#Author Meenakshi Venkatasubramanian - altanalyze@gmail.com

#Permission is hereby granted, free of charge, to any person obtaining a copy 
#of this software and associated documentation files (the "Software"), to deal 
#in the Software without restriction, including without limitation the rights 
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
#copies of the Software, and to permit persons to whom the Software is furnished 
#to do so, subject to the following conditions:

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

""" MetaSpliceWorkflow Module
https://github.com/venkatmi/oncosplice
Steps applied in this workflow:
1 - Run splice-ICGS (Feature Selection)
2 - Block identification (Rank analysis)
3 - NMF Analysis (Initial subtype identification)
4 - Filter Event Annotation
5 - Meta data analysis (differential expression)
6 - Expand clusters (SVM sample classification)
7 - Mutation enrichment (MAF or VCF - optional)
8 - Correlation depletion (excluded biological confounding signatures)
"""
import sys, string, os
import RNASeq
import RNASeq_blockIdentification
import NMF_Analysis as NMF_Analysis
import filterEventAnnotation
import metaDataAnalysis
import ExpandSampleClusters as ExpandSampleClusters
import sampleIndexSelection
import Correlationdepletion
import UI
import multiprocessing as mlp
import export
upd_guides=[]
import operator
from collections import OrderedDict
from collections import defaultdict
import Kmeans
import MutationEnrichment_adj as ME
import Orderedheatmap

def filterPSIValues(filename):
    #fn = filepath(filename)
    firstRow=True
    #detected=0.75   
    header = True
    rows=0
    filtered=0
    new_file = filename[:-4]+'-75p.txt'
   
    ea = export.ExportFile(new_file)

    for line in open(filename,'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data,'\t')
        if header:
            header = False
            eventindex=t.index('EventAnnotation')
            t = [t[1]]+t[eventindex+1:]
            header_length = len(t)-1
            minimum_values_present = int(float((header_length)-1.0)*0.75)
            not_detected = header_length-minimum_values_present
            new_line = line
            ea.write(new_line)
        else:
          
            t = [t[1]]+t[eventindex+1:]
            missing_values_at_the_end = (header_length+1)-len(t)
            missing = missing_values_at_the_end+t.count('')
            if missing<not_detected:
               
                new_line = line
                ea.write(new_line)
               
                filtered+=1
        rows+=1

    ea.close()
    return new_file
    

def header_list(EventAnnot):
    head=0
    header=[]
    with open(EventAnnot, 'rU') as fin:
        for line in fin:
            if head==0:
                line = line.rstrip(os.linesep)
                line=string.split(line,'\t')
                startpos=line.index('EventAnnotation')
                header.append('UID')
                for i in range(startpos+1,len(line)):
                    
                        header.append(line[i])
                
              
                head=1
            else:break
    return header


def FindTopUniqueEvents(Guidefile,psi,Guidedir):
    head=0
    guidekeys=[]
    exportnam=os.path.join(Guidedir,"SplicingeventCount1.txt")
    export_class=open(exportnam,"a")

    tempkeys={}
    global upd_guides
    global train
    omitcluster=0
    
    unique_clusters={}

    for line in open(Guidefile,'rU').xreadlines():
        if head==0:
            line1=line.rstrip('\r\n')
            q= string.split(line1,'\t')
            head=1
            try:
                uid=q.index('UID')
                adjp=q.index('rawp')
                dpsi=q.index('dPSI')
                Clusterid=q.index('UpdatedClusterID')
                cutoff=0.1
                continue
            except Exception:
                uid=q.index('GeneID')
                adjp=q.index('rawp')
                dpsi=q.index('LogFold')
                Clusterid=q.index('GeneID')
                cutoff=0.58
        else:
            line1=line.rstrip('\r\n')
            q= string.split(line1,'\t')
            if abs(float(q[dpsi]))>cutoff and float(q[adjp])<0.01:
                try:
                    tempkeys[q[Clusterid]].append([q[uid],float(q[adjp]),q[adjp+1]])
                except KeyError:
                    tempkeys[q[Clusterid]]=[[q[uid],float(q[adjp]),q[adjp+1]],]
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
    
    try:
        if len(unique_clusters[0])>1:     
            unique_clusters[0].sort(key=operator.itemgetter(1))
            if len(unique_clusters[0])>120:
                guidekeys=unique_clusters[0][0:150]
                for i in range(0,len(guidekeys)):
                    upd_guides.append(guidekeys[i][0])
            else:
                        omitcluster=1
        else:
            omitcluster=1
        export_class.write(psi+"\t"+str(len(unique_clusters[0]))+"\n")
    except Exception:
        omitcluster=1
    
    return omitcluster

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def MergeResults(dire):
    file_index={}
    count=0
    for filename in os.listdir(dire):

        if ("_Results" in filename or "Kmeans" in filename)  and "._" not in filename and "ordered" not in filename:
            file_index[filename]=count
            count+=1

    keylist={}
    
    heads={}
    for filename in os.listdir(dire):
        if ("_Results" in filename or "Kmeans" in filename)  and "._" not in filename and "ordered" not in filename:
            Guidefile=os.path.join(dire, filename)
            head=0
            for line in open(Guidefile,'rU').xreadlines():
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                header=[]
                if head==0:
                    head=1
                    for i in range(1,len(t)):
                        header.append(t[i])
                    heads[filename]=header
                        
                    continue
                else:
                    
                    val=[]
                    key=t[0]
                    for i in range(1,len(t)):
                        val.append(t[i])
                    if key not in keylist:
                        keylist[key]=[[file_index[filename],val],]
                    else:
                        keylist[key].append([file_index[filename],val])
    exportnam=os.path.join(dire,"MergedResult.txt")
    export_class=open(exportnam,"w")
    export_class.write("uid")
    
    for filename in file_index:
        export_class.write("\t")
        print filename,heads[filename]
        export_class.write(string.join(heads[filename],"\t"))
    export_class.write("\n")

    for key in keylist:
        export_class.write(key)
        
        for filename in file_index:
            
            for val1,val2 in keylist[key]:
                if file_index[filename]==val1:
                    export_class.write("\t")
                    export_class.write(string.join(val2,"\t"))
               
                    break
        export_class.write("\n")  
    return exportnam
def CompleteWorkflow(InputFile,EventAnnot,turn,rho_cutoff,strategy,seq):
        
        species="Hs"
        row_method = 'hopach'
        column_method = 'hopach'
        row_metric = 'correlation'
        column_metric = 'euclidean'
        color_gradient = 'yellow_black_blue'
        contrast=3
        vendor = "RNASeq"
        GeneSelection = ''
        PathwaySelection = ''
        GeneSetSelection = 'None Selected'
        excludeCellCycle = False
        #rho_cutoff = 0.4
        restrictBy = 'protein_coding'
        featurestoEvaluate = 'Genes'
        ExpressionCutoff = 0
        CountsCutoff = 0
        FoldDiff = 1.2
        SamplesDiffering = 4
        JustShowTheseIDs=''
        removeOutliers = False
        PathwaySelection=[]
        array_type="RNASeq"
        #rho_cutoff=0.4
        gsp = UI.GeneSelectionParameters(species,array_type,vendor)
        gsp.setGeneSet(GeneSetSelection)
        gsp.setPathwaySelect(PathwaySelection)
        gsp.setGeneSelection(GeneSelection)
        gsp.setJustShowTheseIDs(JustShowTheseIDs)
        gsp.setNormalize('median')
        gsp.setSampleDiscoveryParameters(ExpressionCutoff,CountsCutoff,FoldDiff,SamplesDiffering,removeOutliers,featurestoEvaluate,restrictBy,excludeCellCycle,column_metric,column_method,rho_cutoff) 
        #Run splice ICGS
        
        """import UI
        species='Mm'; platform = "3'array"; vendor = 'Ensembl'
        gsp = UI.GeneSelectionParameters(species,platform,vendor)
        gsp.setGeneSet('None Selected')
        gsp.setPathwaySelect('')
        gsp.setGeneSelection('')
        gsp.setJustShowTheseIDs('')
        gsp.setNormalize('median')
        gsp.setSampleDiscoveryParameters(0,0,1.5,3,
        False,'PSI','protein_coding',False,'cosine','hopach',0.35)"""
        
        FilteredEventAnnot=filterEventAnnotation.FilterFile(InputFile,EventAnnot,turn)
       
        
        try:
            print "Running splice-ICGS for feature selection - Round"+str(turn)
        #except Exception:Rank=0
            graphic_links3 = RNASeq.singleCellRNASeqWorkflow(species, 'exons', InputFile,mlp,exp_threshold=0, rpkm_threshold=0, parameters=gsp)

            Guidefile=graphic_links3[-1][-1]
            Guidefile=Guidefile[:-4]+'.txt'
           
            
            print "Running block identification for rank analyses - Round"+str(turn)
            RNASeq_blockIdentification.correlateClusteredGenesParameters(Guidefile,rho_cutoff=0.4,hits_cutoff=4,hits_to_report=50,ReDefinedClusterBlocks=True,filter=True) 
            Guidefile_block=Guidefile[:-4]+'-BlockIDs.txt'
            NMFinput,Rank=NMF_Analysis.FilterFile(Guidefile,Guidefile_block,InputFile,turn)
        except Exception:Rank=0
     
           

        if Rank>1:
         
            if Rank>2:Rank=30
            else: Rank=2
            if seq=="bulk":
               use_adjusted_p=True
            else:
               use_adjusted_p=False
         
            print "Running NMF analyses for dimension reduction using "+str(Rank)+" ranks - Round"+str(turn)
            NMFResult,BinarizedOutput,Metadata,Annotation=NMF_Analysis.NMFAnalysis(NMFinput,Rank,turn,strategy)
            print "Running Metadata Analyses for finding differential splicing events"
            rootdir,CovariateQuery=metaDataAnalysis.remoteAnalysis('Hs',FilteredEventAnnot,Metadata,'PSI',0.1,use_adjusted_p,0.05,Annotation)
            counter=1
            Guidedir=rootdir+CovariateQuery
            PSIdir=rootdir+'ExpressionProfiles'
            global upd_guides
            upd_guides=[]
            name=[]
            group=[]
            grplst=[]
            for filename in os.listdir(Guidedir):
                if filename.startswith("PSI."):
                    Guidefile=os.path.join(Guidedir, filename)
                    psi=string.replace(filename,"PSI.","")
                    PSIfile=os.path.join(PSIdir, psi)
                    omitcluster=FindTopUniqueEvents(Guidefile,psi,Guidedir)
                   
                    if omitcluster==0:
                        group.append(counter)
                        name.append(psi)
                        counter+=1
            if counter>2:
                dire = export.findParentDir(InputFile)
                output_dir = dire+'OncoInputs'
                if os.path.exists(output_dir)==False:
                    export.createExportFolder(output_dir)
        
                output_file = output_dir+'/SVMInput-Round'+str(turn)+'.txt'
                ExpandSampleClusters.filterRows(InputFile,output_file,filterDB=upd_guides,logData=False)
                header=ExpandSampleClusters.header_file(output_file)
                print "Running SVM prediction for improved subtypes - Round"+str(turn)
                train=ExpandSampleClusters.TrainDataGeneration(output_file,BinarizedOutput,name)
                grplst.append(group)
                ExpandSampleClusters.Classify(header,train,output_file,grplst,name,turn)
                header=Correlationdepletion.header_file(NMFResult)
                
                output_file=output_dir+'/DepletionInput-Round'+str(turn)+".txt"
                sampleIndexSelection.filterFile(InputFile,output_file,header)
                print "Running Correlation Depletion - Round"+str(turn)
                commonkeys,count=Correlationdepletion.FindCorrelations(NMFResult,output_file,name)
                Depleted=Correlationdepletion.DepleteSplicingevents(commonkeys,output_file,count,InputFile)
                InputFile=Depleted
            
                flag=True
            else:
                try:
                    print "Running K-means analyses instead of NMF - Round"+str(turn)
                    header=[]
                    header=Kmeans.header_file(Guidefile_block)
                    Kmeans.KmeansAnalysis(Guidefile_block,header,InputFile,turn)
                    flag=False
                except Exception:
                    flag=False
                
        else:
            if Rank==1:
                try:
                    print "Running K-means analyses instead of NMF - Round"+str(turn)
                    header=[]
                    header=Kmeans.header_file(Guidefile_block)
                    Kmeans.KmeansAnalysis(Guidefile_block,header,InputFile,turn)
            
                    flag=False
                except Exception:
                    flag=False
            else:
                flag=False
         
        return flag,InputFile,FilteredEventAnnot
if __name__ == '__main__':
    import getopt
    seq="bulk"
    rho_cutoff=0.4
    strategy="stringent"
    filters=True
    mode="iterative"
    Mutationref=""
    flag=True
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Insufficient command line flags supplied."
        sys.exit()
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['EventAnnotation=','rho=','strategy=','filter=','mode=','Mutationref=','Assoc='])
        for opt, arg in options:
            #if opt == '--InputFile': InputFile=arg
            if opt=='--EventAnnotation':EventAnnot=arg
            if opt=='--rho':rho_cutoff=arg
            if opt=='--strategy':strategy=arg
            if opt=='--filter':filters=arg
            if opt=='--mode':mode=arg
            if opt=='--Mutationref':Mutationref=arg
            if opt=='--Assoc':assoc=arg
   
    if assoc=="False": flag=False
    
    print EventAnnot
    print Mutationref
    dire = export.findParentDir(EventAnnot)
    turn=1
    if turn==1 and filters==True:
        EventAnnot=filterPSIValues(EventAnnot)
    output_dir = dire+'ExpressionInput'
   
    export.createExportFolder(output_dir)
    InputFile=output_dir+"/exp.input.txt"
    header=header_list(EventAnnot)
    sampleIndexSelection.filterFile(EventAnnot,InputFile,header,FirstCol=False)
    
    if flag: 
        if mode=="single":
            flag,InputFile,EventAnnot=CompleteWorkflow(InputFile,EventAnnot,turn,rho_cutoff,strategy,seq)
      
        else:
            while flag:
                
                flag,InputFile,EventAnnot=CompleteWorkflow(InputFile,EventAnnot,turn,rho_cutoff,strategy,seq)
          
                turn+=1
                if turn>3:
                    flag=False
                if flag==False:
                    break
        
        output_dir = dire+'SVMOutputs'
        Combinedres=MergeResults(output_dir)
    else:
        output_dir = dire+'SVMOutputs'
        Combinedres=os.path.join(output_dir,"MergedResult.txt")
    mutlabels={}
    if Mutationref!="":
        print "Running Mutation Enrichment Analyses"
        Expand="yes"
        mutdict=defaultdict(list)
        header=ME.header_file(Mutationref)
      
        mutdict=ME.findsiggenepermut(Mutationref)
      
        mutlabels=ME.Enrichment(Combinedres,mutdict,Mutationref,Expand,header)
    print "Generating the final consolidated results"
    Orderedheatmap.Classify(Combinedres,mutlabels,dire)
    Orderedheatmap.Classify(Combinedres,mutlabels,dire,False) 
    
    print "successfully completed"
