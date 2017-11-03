#!/usr/bin/env python

import sys, string, os
import RNASeq
import RNASeq_blockIdentification
import NMF_Analysis
import filterEventAnnotation
import metaDataAnalysis
import ExpandClusters_v3_KD as ExpandSampleClusters
import sampleIndexSelection
import Correlationdepletion
import UI
import multiprocessing as mlp
import export
upd_guides=[]
import operator
from collections import OrderedDict
from collections import defaultdict


def FindTopUniqueEvents(Guidefile,psi,Guidedir):
    head=0
    guidekeys=[]
    exportnam=os.path.join(Guidedir,"SplicingeventCount1.txt")
    export_class=open(exportnam,"a")
    #commonkeys=[]
    tempkeys={}
    #upd_guides=[]
    global upd_guides
    train=[]
    omitcluster=0
    
    unique_clusters={}

    for line in open(Guidefile,'rU').xreadlines():
        if head==0:
            head=1
            continue
        else:
            line1=line.rstrip('\r\n')
            q= string.split(line1,'\t')
          
            #if abs(float(q[8]))>0.1 and float(q[10])<0.05:
            if abs(float(q[9]))>0.1:
                try:
                    tempkeys[q[3]].append([q[0],float(q[10]),q[12]])
                except KeyError:
                    tempkeys[q[3]]=[[q[0],float(q[10]),q[12]],]
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
  
            if len(unique_clusters[0])>10:
                guidekeys=unique_clusters[0][0:300]
                for i in range(0,len(guidekeys)):
            
            #upd_guides[i]=[upd_guides[i][3],upd_guides[i][4]]
                    upd_guides.append(guidekeys[i][0])
            else:
            #        #if len(unique_clusters[0])>50:
            #            guidekeys=unique_clusters[0][0:25]
            #            for i in range(0,len(guidekeys)):
            #
            ##upd_guides[i]=[upd_guides[i][3],upd_guides[i][4]]
            #                upd_guides.append(guidekeys[i][0])
                    #else:
                        omitcluster=1
            export_class.write(psi+"\t"+str(len(unique_clusters[0]))+"\n")
    except Exception:
        omitcluster=1
    
    print len(upd_guides)
    return omitcluster
def CompleteWorkflow(InputFile,EventAnnot,turn):
    
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
        rho_cutoff = 0.4
        restrictBy = 'protein_coding'
        featurestoEvaluate = 'AltExon'
        ExpressionCutoff = 0
        CountsCutoff = 0
        FoldDiff = 1.3
        SamplesDiffering = 4
        JustShowTheseIDs=''
        removeOutliers = False
        PathwaySelection=[]
        array_type="RNASeq"
        rho_cutoff=0.4
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
        
        Rank=30
        if Rank>1:
            if Rank>3:Rank=30
            else: Rank=2
       
            FilteredEventAnnot=EventAnnot
   
            counter=1
            #Guidedir=rootdir+CovariateQuery
            #PSIdir=rootdir+'ExpressionProfiles'
       
            #BinarizedOutput="/Volumes/Pass/KnownClasses/AllAnnotations_T.txt"
            Guidedir="/Volumes/Pass/ConsolidatedResults/MLL2"
            PSIdir="/Volumes/Pass/MOdifiedNMF/round1_old/ExpressionProfiles"
            Guidedir="/Volumes/Pass/Knockdown_Classification/RBP-KD/ENCODE/Events-dPSI_0.1_rawp"
            PSIdir="/Volumes/Pass/Knockdown_Classification/RBP-KD/ENCODE/ExpressionProfiles"
            BinarizedOutput="/Volumes/Pass/Knockdown_Classification/RBP-KD/ENCODE/ExpressionInput/Binarized_Output.txt"
            #BinarizedOutput="/Volumes/Pass/MOdifiedNMF/ExpressionInput/exp.splicing-filteredcor_depleted-filteredNMFsnmf_binary30_filtered.txt"
            global upd_guides
            upd_guides=[]
            name=[]
            group=[]
            grplst=[]
            filteredevents=[]
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
            #InputFile2=""
            #InputFile2="/Volumes/Pass/Leucegene_ClassifyNew/exp.splicing-filteredcor_depleted.txt"
            InputFile2="/Volumes/Pass/Knockdown_Classification/RBP-KD/ENCODE/ExpressionInput/exp.RBP-KDs.txt"
            output_file=InputFile[:-4]+"-filtered.txt"
            output_file1=InputFile2[:-4]+"-filtered.txt" 
            #print guidekey
            #print upd_guides
            print len(upd_guides)
            filteredevents=ExpandSampleClusters.filterRows_data(InputFile,output_file,filterDB=upd_guides,logData=False)
            ExpandSampleClusters.filterRows_data(InputFile2,output_file1,filterDB=filteredevents,logData=False)
            header=ExpandSampleClusters.header_file(output_file)
            train,train2=ExpandSampleClusters.TrainDataGeneration(output_file1,BinarizedOutput,name)
            grplst.append(group)
            Finalclusters,FinalCentroids=ExpandSampleClusters.Classify(header,train,train2,output_file,grplst,name)
            #header=Correlationdepletion.header_file(NMFResult)
            #output_file=InputFile[:-4]+"-filtered.txt"
            #sampleIndexSelection.filterFile(InputFile,output_file,header)
            #commonkeys,count=Correlationdepletion.FindCorrelations(NMFResult,output_file)
            #Depleted=Correlationdepletion.DepleteSplicingevents(commonkeys,output_file,count)
            #InputFile=Depleted
            
        #    flag=True
        #else:
        #    flag=False
        ##    
        #    
        #return flag,InputFile,FilteredEventAnnot
if __name__ == '__main__':
    import getopt
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Insufficient command line flags supplied."
        sys.exit()
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['InputFile=','EventAnnotation='])
        for opt, arg in options:
            if opt == '--InputFile': InputFile=arg
            elif opt=='--EventAnnotation':EventAnnot=arg
    flag=True
    turn=1
    while flag:
        flag,InputFile,EventAnnot=CompleteWorkflow(InputFile,EventAnnot,turn)
        turn+=1
        flag=False
        if flag==False:
            break
    print "completed"