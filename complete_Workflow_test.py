#!/usr/bin/env python


""" Steps involved:
run splice icgs
block identification
NMF Analysis
filter Event Annotation
Meta data analysis
expand clusters
mutation enrichment
correlation depletion
"""
import sys, string, os
import RNASeq
import RNASeq_blockIdentification
import NMF_Analysis_test
import filterEventAnnotation
import metaDataAnalysis
import ExpandClusters_v2_test as ExpandSampleClusters
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
            if abs(float(q[9]))>0.1 and float(q[11])<0.05:
                try:
                    
                    tempkeys[q[3]].append([q[0],float(q[11]),q[12]])
                except KeyError:
                    tempkeys[q[3]]=[[q[0],float(q[11]),q[12]],]
            
    for i in tempkeys:
        
        if len(tempkeys[i])>1:
          for ijq in range(0,len(tempkeys[i])):
        
          
            #print tempkeys[i]
            tempkeys[i].sort(key=operator.itemgetter(1),reverse=False)
            #print tempkeys[i][0]
            try:
                #unique_clusters[0].append(tempkeys[i])
               unique_clusters[0].append(tempkeys[i][ijq])
            except KeyError:
                #unique_clusters[0]=[tempkeys[i],]
                unique_clusters[0]=[tempkeys[i][ijq],]
        else:
            try:
                unique_clusters[0].append(tempkeys[i][0])
            except KeyError:
                unique_clusters[0]=[tempkeys[i][0],]
    
    try:
        print len(unique_clusters[0])
        if len(unique_clusters[0])>1:     
            unique_clusters[0].sort(key=operator.itemgetter(1))
  
            #if len(unique_clusters[0])>150:
            guidekeys=unique_clusters[0][0:30]
            
              #  for i in range(0,len(guidekeys)):
            
            #upd_guides[i]=[upd_guides[i][3],upd_guides[i][4]]
               #     upd_guides.append(guidekeys[i][0])
            #else:
                
            #if len(unique_clusters[0])>99:
            #guidekeys=unique_clusters[0]
            for i in range(0,len(guidekeys)):
            #
            ##upd_guides[i]=[upd_guides[i][3],upd_guides[i][4]]
                upd_guides.append(guidekeys[i][0])
            ##            for i in range(0,len(guidekeys)):
            ##
            ###upd_guides[i]=[upd_guides[i][3],upd_guides[i][4]]
            #                upd_guides.append(guidekeys[i][0])
        else:
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
    
        #Guidefile="/Volumes/Pass/128Samples/ExpressionInput/amplify/DataPlots/Clustering-exp.splicing-filtered-Guide3 LILRA6 ENSG00000244482 E3.3-E4.1 ENSG000002-hierarchical_euclidean_correlation.txt"
    
        #graphic_links3 = RNASeq.singleCellRNASeqWorkflow(species, 'exons', InputFile,mlp,exp_threshold=0, rpkm_threshold=0, parameters=gsp)
        #Guidefile=graphic_links3[-1][-1]
        ##print Guidefile
        #Guidefile=Guidefile[:-4]+'.txt'
        ##Guidefile="/Volumes/Pass/Leucegene_Complete/ICGS/Clustering-exp.splicing-filteredcor_depleted-Guide3 DYNLL1 ENSG00000088986 E2.1-I2.1 ENSG000000-hierarchical_euclidean_correlation.txt"
        #FilteredEventAnnot=filterEventAnnotation.FilterFile(InputFile,EventAnnot,turn)
        ##Guidefile='/Users/meenakshi/Documents/testdata/ExpressionInput/amplify//DataPlots/Clustering-exp.testdata-Guide3 TTC27 ENSG00000018699 I1.1_32854940-E2.1 EN-hierarchical_euclidean_correlation-BlockIDs.txt'
        #try:
         #   RNASeq_blockIdentification.correlateClusteredGenesParameters(Guidefile,rho_cutoff=0.4,hits_cutoff=4,hits_to_report=50,ReDefinedClusterBlocks=True,filter=True)
        ##    
          #  Guidefile_block=Guidefile[:-4]+'-BlockIDs.txt'
        ###Guidefile="/Users/meenakshi/Documents/testdata/ICGS/Clustering-exp.testdata-filteredcor_depleted-Guide3 ANKS1A ENSG00000064999 I17.1-E18.1 ENSG0000-hierarchical_euclidean_correlation.txt"
        #Guidefile="/Volumes/Pass/NMFTest/Clustering-exp.splicing-filteredcor_depleted-filtered-filteredcor_depleted-Guide3 C16orf62 ENSG00000103544 I15.1-E16.1 ENSG00-hierarchical_euclidean_correlation.txt"
        #Guidefile_block=Guidefile[:-4]+'-BlockIDs.txt'
        #NMFinput,Rank=NMF_Analysis_test.FilterFile(Guidefile,Guidefile_block,InputFile)
        ##    print Rank
        #except Exception:Rank=1
        #sys.exit()
        Rank=30
        if Rank>0:
            if Rank>3:Rank=30
            else: Rank=2
           # NMFResult,BinarizedOutput,Metadata,Annotation=NMF_Analysis_test.NMFAnalysis(NMFinput,Rank,turn)
            #FilteredEventAnnot=filterEventAnnotation.FilterFile(InputFile,EventAnnot,turn)
            #FilteredEventAnnot="/Volumes/Pass/ConsolidatedResults/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt-filtered.txt-filtered.txt"
            #splicingEventTypes={}
            #all_groups_db={}
            #all_comps_db={}
            #platform='PSI'
            #PercentExp=75
            #use_adjusted_p=True
            #CovariateQuery = 'Events'
            #Annotation="/Volumes/Pass/ConsolidatedResults/Annotation.txt"
            #Metadata="/Volumes/Pass/ConsolidatedResults/Allunique_Metadata1.txt"
            #metadata_filters = metaDataAnalysis.importMetaDataDescriptions(Annotation)
            #all_groups_db, all_comps_db = metaDataAnalysis.prepareComparisonData(Metadata,metadata_filters,all_groups_db, all_comps_db)
            ##
            #for i in all_groups_db:
            #    print i
            #    for k in all_groups_db[i]: print '  ',k,'\t',len(all_groups_db[i][k])
            #print all_comps_db
            #
            #if platform == 'PSI':
            #    result_type = 'dPSI'
            #else:
            #    result_type = 'LogFold'
            #logfold_threshold=0.1
            #if use_adjusted_p:
            #    CovariateQuery += '-'+result_type+'_'+str(logfold_threshold)[:4]+'_adjp'
            #else:
            #    CovariateQuery += '-'+result_type+'_'+str(logfold_threshold)+'_rawp'
            #
            #for specificCovariate in all_groups_db:
            #    comps_db = all_comps_db[specificCovariate]
            #    groups_db = all_groups_db[specificCovariate]
            #    rootdir,splicingEventTypes = metaDataAnalysis.performDifferentialExpressionAnalysis(species,platform,FilteredEventAnnot,groups_db,comps_db,CovariateQuery,splicingEventTypes)
            #
            #if platform == 'PSI':
            #    metaDataAnalysis.outputSplicingSummaries(rootdir,splicingEventTypes)
            counter=1
            #Guidedir=rootdir+CovariateQuery
            #PSIdir=rootdir+'ExpressionProfiles'
            ##Guidedir="/Volumes/Pass/MOdifiedNMF/round1/Events-dPSI_0.1_adjp"
            ##PSIdir="/Volumes/Pass/MOdifiedNMF/round1/ExpressionProfiles"
            Guidedir="/Volumes/Pass/MetaData"
            BinarizedOutput='/Volumes/Pass/MetaData/Meta.txt'
            #PSIdir="/Volumes/Pass/Target_SplicingData/round1/ExpressionProfiles"
            #BinarizedOutput="/Volumes/Pass/Target_SplicingData/exp.splicing_filter-filtered-filteredNMFsnmf_binary_specific2.txt"
            global upd_guides
            upd_guides=[]
            name=[]
            group=[]
            grplst=[]
            for filename in os.listdir(Guidedir):
                if filename.startswith("PSI."):
                    Guidefile=os.path.join(Guidedir, filename)
                    psi=string.replace(filename,"PSI.","")
                    #PSIfile=os.path.join(PSIdir, psi)
                    print Guidefile
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
        #    grplst.append(group)
        #    Finalclusters,FinalCentroids=ExpandSampleClusters.Classify(header,train,output_file,grplst,name)
        #    NMFResult=""
        #    header=Correlationdepletion.header_file(NMFResult)
        #    output_file=InputFile[:-4]+"-filtered.txt"
        #    sampleIndexSelection.filterFile(InputFile,output_file,header)
        #    commonkeys,count=Correlationdepletion.FindCorrelations(NMFResult,output_file)
        #    Depleted=Correlationdepletion.DepleteSplicingevents(commonkeys,output_file,count)
        #    InputFile=Depleted
        #    
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
        if flag==False:
            break
    print "completed"