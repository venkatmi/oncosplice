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

def filterPSIValues(filename):
    fn = filepath(filename)
    firstRow=True
          
    header = True
    rows=0
    filtered=0
    new_file = filename[:-4]+'-75p.txt'
    #new_file_clust = new_file[:-4]+'-clustID.txt'
    ea = export.ExportFile(new_file)
    #eac = export.ExportFile(new_file_clust)
   # added=[]
    for line in open(fn,'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data,'\t')
        if header:
            header = False
            eventindex=t.index('EventAnnotation')
            t = [t[1]]+t[eventindex+1:]
            header_length = len(t)-1
            minimum_values_present = int(header_length)-1
            not_detected = header_length-minimum_values_present
            new_line = string.join(t,'\t')+'\n'
            ea.write(new_line)
        else:
            #cID = t[5]
            t = [t[1]]+t[eventindex+1:]
            missing_values_at_the_end = (header_length+1)-len(t)
            missing = missing_values_at_the_end+t.count('')
            if missing<not_detected:
                #if cID not in added:
                #added.append(cID)
                new_line = string.join(t,'\t')+'\n'
                ea.write(new_line)
                #eac.write(t[0]+'\t'+cID+'\n')
                filtered+=1
        rows+=1
    #print rows, filtered
    ea.close()
    return newfile
    #eac.close()
    #removeRedundantCluster(new_file,new_file_clust)

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
                
                #del header[:1]
                head=1
            else:break
    return header


def FindTopUniqueEvents(Guidefile,psi,Guidedir):
    head=0
    guidekeys=[]
    exportnam=os.path.join(Guidedir,"SplicingeventCount1.txt")
    export_class=open(exportnam,"a")
    #commonkeys=[]
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
        if len(unique_clusters[0])>1:     
            unique_clusters[0].sort(key=operator.itemgetter(1))
            if len(unique_clusters[0])>90:
                guidekeys=unique_clusters[0][0:150]
                for i in range(0,len(guidekeys)):
                    upd_guides.append(guidekeys[i][0])
            else:
                        omitcluster=1
            export_class.write(psi+"\t"+str(len(unique_clusters[0]))+"\n")
    except Exception:
        omitcluster=1
    
    print len(upd_guides)
    return omitcluster
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
        
        FilteredEventAnnot=""        
        #try:
        graphic_links3 = RNASeq.singleCellRNASeqWorkflow(species, 'exons', InputFile,mlp,exp_threshold=0, rpkm_threshold=0, parameters=gsp)
        Guidefile=graphic_links3[-1][-1]
        try:
            Guidefile=Guidefile[:-4]+'.txt'
            RNASeq_blockIdentification.correlateClusteredGenesParameters(Guidefile,rho_cutoff=rho_cutoff,hits_cutoff=4,hits_to_report=50,ReDefinedClusterBlocks=True,filter=True)
        ##    
            Guidefile_block=Guidefile[:-4]+'-BlockIDs.txt'
       
        ###Guidefile="/Users/meenakshi/Documents/testdata/ICGS/Clustering-exp.testdata-filteredcor_depleted-Guide3 ANKS1A ENSG00000064999 I17.1-E18.1 ENSG0000-hierarchical_euclidean_correlation.txt"
            NMFinput,Rank=NMF_Analysis.FilterFile(Guidefile,Guidefile_block,InputFile)
           
            FilteredEventAnnot=filterEventAnnotation.FilterFile(InputFile,EventAnnot,turn)
        except Exception:Rank=0
        
        if Rank>1:
         
            if Rank>2:Rank=30
            else: Rank=2
            if seq=="bulk":
               use_adjusted_p=True
            else:
               use_adjusted_p=False
         
           
            NMFResult,BinarizedOutput,Metadata,Annotation=NMF_Analysis.NMFAnalysis(NMFinput,Rank,turn,strategy)
          
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
                    print Guidefile,PSIfile
                    omitcluster=FindTopUniqueEvents(Guidefile,psi,Guidedir)
                    print omitcluster
                    if omitcluster==0:
                        group.append(counter)
                        name.append(psi)
                        counter+=1
        
            output_file=InputFile[:-4]+"-filtered.txt"  
            #print guidekey
            print len(upd_guides)
            ExpandSampleClusters.filterRows(InputFile,output_file,filterDB=upd_guides,logData=False)
            header=ExpandSampleClusters.header_file(output_file)
            train=ExpandSampleClusters.TrainDataGeneration(output_file,BinarizedOutput,name)
            grplst.append(group)
            Finalclusters,FinalCentroids=ExpandSampleClusters.Classify(header,train,output_file,grplst,name)
            header=Correlationdepletion.header_file(NMFResult)
            output_file=InputFile[:-4]+"-filtered.txt"
            sampleIndexSelection.filterFile(InputFile,output_file,header)
            commonkeys,count=Correlationdepletion.FindCorrelations(NMFResult,output_file,name)
            Depleted=Correlationdepletion.DepleteSplicingevents(commonkeys,output_file,count)
            InputFile=Depleted
            
            flag=True
        else:
            try:
                header=[]
                header=Kmeans.header_file(Guidefile_block)
                Kmeans.KmeansAnalysis(InputFile,Guidefile_block,header)
            
                flag=False
            except Exception:
                flag=False
        #    
            
        return flag,InputFile,FilteredEventAnnot
if __name__ == '__main__':
    import getopt
    seq="bulk"
    rho_cutoff=0.4
    strategy="stringent"
    filters=False
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Insufficient command line flags supplied."
        sys.exit()
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['EventAnnotation=','rho=','strategy=','seq=','filter='])
        for opt, arg in options:
            #if opt == '--InputFile': InputFile=arg
            if opt=='--EventAnnotation':EventAnnot=arg
            if opt=='--rho':rho_cutoff=arg
            if opt=='--strategy':strategy=arg
            if opt=='--seq':seq=arg
            if opt=='--filter':filters=arg
            
    dire = export.findParentDir(EventAnnot)
   
    output_dir = dire+'ExpressionInput'
    export.createExportFolder(output_dir)
    InputFile=output_dir+"/exp.input.txt"
    
    header=header_list(EventAnnot)
   
    sampleIndexSelection.filterFile(EventAnnot,InputFile,header,FirstCol=False)
    flag=True
    turn=1
    
    while flag:
        if turn==1 and filters==True:
            InputFile=filterPSIValues(InputFile)
       
        flag,InputFile,EventAnnot=CompleteWorkflow(InputFile,EventAnnot,turn,rho_cutoff,strategy,seq)
      
        turn+=1
        if flag==False:
            break
    print "completed"