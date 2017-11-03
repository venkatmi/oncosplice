###RNASeq
#Copyright 2005-2008 J. David Gladstone Institutes, San Francisco California
#Author Nathan Salomonis - nsalomonis@gmail.com

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

import sys, string
import statistics
import math
import os.path
import unique
import update
import copy
import time
import export
import EnsemblImport; reload(EnsemblImport)
import JunctionArrayEnsemblRules
import JunctionArray; reload(JunctionArray)
import ExonArrayEnsemblRules
import multiprocessing
import logging
import traceback
import warnings

try:
    from scipy import average as Average
except Exception:
    from statistics import avg as Average
    
def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    return dir_list

def makeUnique(item):
    db1={}; list1=[]; k=0
    for i in item:
        try: db1[i]=[]
        except TypeError: db1[tuple(i)]=[]; k=1
    for i in db1:
        if k==0: list1.append(i)
        else: list1.append(list(i))
    list1.sort()
    return list1

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

######### Below code deals with building the AltDatabase #########
def collapseNoveExonBoundaries(novel_exon_coordinates,dataset_dir):
    """ Merge exon predictions based on junction measurments from TopHat. The predicted exons are
    bound by the identified splice site and the consensus length of reads in that sample"""
    
    dataset_dir = string.replace(dataset_dir,'exp.','ExpressionInput/novel.')
    
    export_data,status = AppendOrWrite(dataset_dir) ### Export all novel exons
    if status == 'not found':
        export_data.write('GeneID\tStrand\tExonID\tCoordinates\n')
    novel_gene_exon_db={}
    for (chr,coord) in novel_exon_coordinates:
        key = (chr,coord)
        ji,side,coord2 = novel_exon_coordinates[(chr,coord)]
        try:
            if side == 'left': ### left corresponds to the position of coord
                intron = string.split(ji.ExonRegionID(),'-')[1][:2]
            else:
                intron = string.split(ji.ExonRegionID(),'-')[0][:2]
            ls = [coord,coord2]
            ls.sort() ### The order of this is variable
            if ji.Strand() == '-':
                coord2,coord = ls
            else: coord,coord2 = ls
            if 'I' in intron and ji.Novel() == 'side':
                #if 'ENSG00000221983' == ji.GeneID():
                try: novel_gene_exon_db[ji.GeneID(),ji.Strand(),intron].append((coord,coord2,ji,key,side))
                except Exception: novel_gene_exon_db[ji.GeneID(),ji.Strand(),intron] = [(coord,coord2,ji,key,side)]
        except Exception: pass
    
    outdatedExons={} ### merging novel exons, delete one of the two original
    for key in novel_gene_exon_db:
        firstNovel=True ### First putative novel exon coordinates examined for that gene
        novel_gene_exon_db[key].sort()
        if key[1]=='-':
            novel_gene_exon_db[key].reverse()

        for (c1,c2,ji,k,s) in novel_gene_exon_db[key]:
            if firstNovel==False:
                #print [c1,l2] #abs(c1-l2);sys.exit()
                ### see if the difference between the start position of the second exon is less than 300 nt away from the end of the last
                if abs(c2-l1) < 300 and os!=s: ### 80% of human exons are less than 200nt - PMID: 15217358
                    proceed = True
                    #if key[1]=='-':
                    if c2 in k:
                        novel_exon_coordinates[k] = ji,s,l1
                        outdatedExons[ok]=None ### merged out entry
                    elif l1 in ok:
                        novel_exon_coordinates[ok] = li,os,c2
                        outdatedExons[k]=None ### merged out entry
                    else:
                        proceed = False ### Hence, the two splice-site ends are pointing to two distinct versus one common exons
                    """
                    if c2 == 18683670 or l1 == 18683670:
                        print key,abs(c2-l1), c1, c2, l1, l2, li.ExonRegionID(), ji.ExonRegionID();
                        print k,novel_exon_coordinates[k]
                        print ok,novel_exon_coordinates[ok]
                    """
                    if proceed:
                        values = string.join([ji.GeneID(),ji.Strand(),key[2],ji.Chr()+':'+str(l1)+'-'+str(c2)],'\t')+'\n'
                        export_data.write(values)
                    
            ### For negative strand genes, c1 is larger than c2 but is the 5' begining of the exon
            l1,l2,li,ok,os = c1,c2,ji,k,s ### record the last entry
            firstNovel=False
            
    for key in outdatedExons: ### Delete the non-merged entry
        del novel_exon_coordinates[key]
    export_data.close()
    return novel_exon_coordinates
    

def exportNovelExonToBedCoordinates(species,novel_exon_coordinates,chr_status,searchChr=None):
    ### Export the novel exon coordinates based on those in the junction BED file to examine the differential expression of the predicted novel exon
    #bamToBed -i accepted_hits.bam -split| coverageBed -a stdin -b /home/databases/hESC_differentiation_exons.bed > day20_7B__exons-novel.bed
    bed_export_path = filepath('AltDatabase/'+species+'/RNASeq/chr/'+species + '_Ensembl_exons'+searchChr+'.bed')
    bed_data = open(bed_export_path,'w') ### Appends to existing file
    for (chr,coord) in novel_exon_coordinates:
        ji,side,coord2 = novel_exon_coordinates[(chr,coord)]
        if side == 'left': start,stop = coord,coord2
        if side == 'right': start,stop = coord2,coord
        try: gene = ji.GeneID()
        except Exception: gene = 'NA'
        if gene == None: gene = 'NA'
        if gene == None: gene = 'NA'
        if gene != 'NA': ### Including these has no benefit for AltAnalyze (just slows down alignment and piles up memory)
            if ji.Strand() == '-': stop,start=start,stop
            if chr_status == False:
                chr = string.replace(chr,'chr','') ### This will thus match up to the BAM files
            bed_values = [chr,str(start),str(stop),gene,'0',str(ji.Strand())]
            bed_values = cleanUpLine(string.join(bed_values,'\t'))+'\n'
            bed_data.write(bed_values)
    bed_data.close()
    return bed_export_path

def moveBAMtoBEDFile(species,dataset_name,root_dir):
    bed_export_path = filepath('AltDatabase/'+species+'/RNASeq/'+species + '_Ensembl_exons.bed')
    dataset_name = string.replace(dataset_name,'exp.','')
    new_fn = root_dir+'/BAMtoBED/'+species + '_'+dataset_name+'_exons.bed'
    new_fn = string.replace(new_fn,'.txt','')
    print 'Writing exon-level coordinates to BED file:'
    print new_fn
    catFiles(bed_export_path,'chr') ### concatenate the files ot the main AltDatabase directory then move
    export.customFileMove(bed_export_path,new_fn)
    return new_fn

def reformatExonFile(species,type,chr_status):
    if type == 'exon':
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_exon.txt'
        export_path = 'AltDatabase/'+species+'/RNASeq/'+species + '_Ensembl_exons.txt'
        ### Used by BEDTools to get counts per specific AltAnalyze exon region (should augment with de novo regions identified from junction analyses)
        bed_export_path = 'AltDatabase/'+species+'/RNASeq/chr/'+species + '_Ensembl_exons.bed'
        bed_data = export.ExportFile(bed_export_path)
    else:
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_junction.txt'
        export_path = 'AltDatabase/'+species+'/RNASeq/'+species + '_Ensembl_junctions.txt'
    print 'Writing',export_path
    export_data = export.ExportFile(export_path)
    fn=filepath(filename); x=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0:
            x+=1
            export_title = ['AltAnalyzeID','exon_id','ensembl_gene_id','transcript_cluster_id','chromosome','strand','probeset_start','probeset_stop']
            export_title +=['affy_class','constitutive_probeset','ens_exon_ids','ens_constitutive_status','exon_region','exon-region-start(s)','exon-region-stop(s)','splice_events','splice_junctions']
            export_title = string.join(export_title,'\t')+'\n'; export_data.write(export_title)
        else:
            gene, exonid, chr, strand, start, stop, constitutive_call, ens_exon_ids, splice_events, splice_junctions = t
            if chr == 'chrM': chr = 'chrMT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            if chr == 'M': chr = 'MT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention,
            if constitutive_call == 'yes': ens_constitutive_status = '1'
            else: ens_constitutive_status = '0'
            export_values = [gene+':'+exonid, exonid, gene, '', chr, strand, start, stop, 'known', constitutive_call, ens_exon_ids, ens_constitutive_status]
            export_values+= [exonid, start, stop, splice_events, splice_junctions]
            export_values = string.join(export_values,'\t')+'\n'; export_data.write(export_values)
            if type == 'exon':
                if chr_status == False:
                    chr = string.replace(chr,'chr','') ### This will thus match up to the BAM files
                bed_values = [chr,start,stop,gene+':'+exonid+'_'+ens_exon_ids,'0',strand]
                bed_values = string.join(bed_values,'\t')+'\n'; bed_data.write(bed_values)
    export_data.close()
    if type == 'exon': bed_data.close()

def importExonAnnotations(species,type,search_chr):
    if 'exon' in type:
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_exon.txt'
    else:
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_Ensembl_junction.txt'

    fn=filepath(filename); x=0; exon_annotation_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: x=1
        else:
            gene, exonid, chr, strand, start, stop, constitutive_call, ens_exon_ids, splice_events, splice_junctions = t; proceed = 'yes'
            if chr == 'chrM': chr = 'chrMT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            if chr == 'M': chr = 'MT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            if len(search_chr)>0:
                if chr != search_chr: proceed = 'no'
            if proceed == 'yes':
                if type == 'exon': start = int(start); stop = int(stop)
                ea = EnsemblImport.ExonAnnotationsSimple(chr, strand, start, stop, gene, ens_exon_ids, constitutive_call, exonid, splice_events, splice_junctions)
                if type == 'junction_coordinates': 
                    exon1_start,exon1_stop = string.split(start,'|')
                    exon2_start,exon2_stop = string.split(stop,'|')
                    if strand == '-':
                        exon1_stop,exon1_start = exon1_start,exon1_stop
                        exon2_stop,exon2_start = exon2_start,exon2_stop
                    #if gene == 'ENSMUSG00000027340': print chr,int(exon1_stop),int(exon2_start)
                    exon_annotation_db[chr,int(exon1_stop),int(exon2_start)]=ea              
                elif type == 'distal-exon':
                    exon_annotation_db[gene] = exonid
                else:
                    try: exon_annotation_db[gene].append(ea)
                    except KeyError: exon_annotation_db[gene]=[ea]
    return exon_annotation_db

def exportKnownJunctionComparisons(species):
    gene_junction_db = JunctionArrayEnsemblRules.importEnsemblUCSCAltJunctions(species,'standard')
    gene_intronjunction_db = JunctionArrayEnsemblRules.importEnsemblUCSCAltJunctions(species,'_intronic')
    for i in gene_intronjunction_db: gene_junction_db[i]=[]

    gene_junction_db2={}
    for (gene,critical_exon,incl_junction,excl_junction) in gene_junction_db:
        critical_exons = string.split(critical_exon,'|')
        for critical_exon in critical_exons:
            try: gene_junction_db2[gene,incl_junction,excl_junction].append(critical_exon)
            except Exception: gene_junction_db2[gene,incl_junction,excl_junction] = [critical_exon]
    gene_junction_db = gene_junction_db2; gene_junction_db2=[]
    
    junction_export = 'AltDatabase/' + species + '/RNASeq/'+ species + '_junction_comps.txt'
    fn=filepath(junction_export); data = open(fn,'w')
    print "Exporting",junction_export
    title = 'gene'+'\t'+'critical_exon'+'\t'+'exclusion_junction_region'+'\t'+'inclusion_junction_region'+'\t'+'exclusion_probeset'+'\t'+'inclusion_probeset'+'\t'+'data_source'+'\n'
    data.write(title); temp_list=[]
    for (gene,incl_junction,excl_junction) in gene_junction_db:
        critical_exons = unique.unique(gene_junction_db[(gene,incl_junction,excl_junction)])
        critical_exon = string.join(critical_exons,'|')
        temp_list.append(string.join([gene,critical_exon,excl_junction,incl_junction,gene+':'+excl_junction,gene+':'+incl_junction,'AltAnalyze'],'\t')+'\n')
    temp_list = unique.unique(temp_list)
    for i in temp_list: data.write(i)
    data.close()

def getExonAndJunctionSequences(species):
    export_exon_filename = 'AltDatabase/'+species+'/RNASeq/'+species+'_Ensembl_exons.txt'    
    ensembl_exon_db = ExonArrayEnsemblRules.reimportEnsemblProbesetsForSeqExtraction(export_exon_filename,'null',{})

    ### Import just the probeset region for mRNA alignment analysis
    analysis_type = ('region_only','get_sequence'); array_type = 'RNASeq'
    dir = 'AltDatabase/'+species+'/SequenceData/chr/'+species; gene_seq_filename = dir+'_gene-seq-2000_flank.fa'
    ensembl_exon_db = EnsemblImport.import_sequence_data(gene_seq_filename,ensembl_exon_db,species,analysis_type)

    critical_exon_file = 'AltDatabase/'+species+'/'+ array_type + '/' + array_type+'_critical-exon-seq.txt'
    getCriticalJunctionSequences(critical_exon_file,species,ensembl_exon_db)

    """    
    ### Import the full Ensembl exon sequence (not just the probeset region) for miRNA binding site analysis   
    analysis_type = 'get_sequence'; array_type = 'RNASeq'
    dir = 'AltDatabase/'+species+'/SequenceData/chr/'+species; gene_seq_filename = dir+'_gene-seq-2000_flank.fa'
    ensembl_exon_db = EnsemblImport.import_sequence_data(gene_seq_filename,ensembl_exon_db,species,analysis_type)
    """
    
    critical_exon_file = 'AltDatabase/'+species+'/'+ array_type + '/' + array_type+'_critical-exon-seq.txt'
    updateCriticalExonSequences(critical_exon_file, ensembl_exon_db)

def updateCriticalExonSequences(filename,ensembl_exon_db):
    exon_seq_db_filename = filename[:-4]+'_updated.txt'
    exonseq_data = export.ExportFile(exon_seq_db_filename)
    
    critical_exon_seq_db={}; null_count={}
    for gene in ensembl_exon_db:
        gene_exon_data={}
        for probe_data in ensembl_exon_db[gene]:
            exon_id,((probe_start,probe_stop,probeset_id,exon_class,transcript_clust),ed) = probe_data
            try: gene_exon_data[probeset_id] = ed.ExonSeq()
            except Exception: null_count[gene]=[] ### Occurs for non-chromosomal DNA (could also download this sequence though)
        if len(gene_exon_data)>0: critical_exon_seq_db[gene] = gene_exon_data
    print len(null_count),'genes not assigned sequenced (e.g.,non-chromosomal)'
    ensembl_exon_db=[]
    
    ### Export exon sequences 
    for gene in critical_exon_seq_db:
        gene_exon_data = critical_exon_seq_db[gene]
        for probeset in gene_exon_data:
            critical_exon_seq = gene_exon_data[probeset]
            values = [probeset,'',critical_exon_seq]
            values = string.join(values,'\t')+'\n'
            exonseq_data.write(values)  
    exonseq_data.close()
    print exon_seq_db_filename, 'exported....'

def getCriticalJunctionSequences(filename,species,ensembl_exon_db):
    ### Assemble and export junction sequences
    junction_seq_db_filename = string.replace(filename,'exon-seq','junction-seq')
    junctionseq_data = export.ExportFile(junction_seq_db_filename)

    critical_exon_seq_db={}; null_count={}
    for gene in ensembl_exon_db:
        gene_exon_data={}
        for probe_data in ensembl_exon_db[gene]:
            exon_id,((probe_start,probe_stop,probeset_id,exon_class,transcript_clust),ed) = probe_data
            try: gene_exon_data[probeset_id] = ed.ExonSeq()
            except Exception: null_count[gene]=[] ### Occurs for non-chromosomal DNA (could also download this sequence though)
        if len(gene_exon_data)>0: critical_exon_seq_db[gene] = gene_exon_data
    print len(null_count),'genes not assigned sequenced (e.g.,non-chromosomal)'
    ensembl_exon_db=[]
    
    junction_annotation_db = importExonAnnotations(species,'junction',[])
    for gene in junction_annotation_db:
        if gene in critical_exon_seq_db:
            gene_exon_data = critical_exon_seq_db[gene]
            for jd in junction_annotation_db[gene]:
                exon1,exon2=string.split(jd.ExonRegionIDs(),'-')
                p1=gene+':'+exon1
                p2=gene+':'+exon2
                p1_seq=gene_exon_data[p1][-15:]
                p2_seq=gene_exon_data[p2][:15]
                junction_seq = p1_seq+'|'+p2_seq
                junctionseq_data.write(gene+':'+jd.ExonRegionIDs()+'\t'+junction_seq+'\t\n')
    junctionseq_data.close()
    print junction_seq_db_filename, 'exported....'

def getEnsemblAssociations(species,data_type,test_status,force):
    ### Get UCSC associations (download databases if necessary)
    import UCSCImport
    mRNA_Type = 'mrna'; run_from_scratch = 'yes'
    export_all_associations = 'no' ### YES only for protein prediction analysis
    update.buildUCSCAnnoationFiles(species,mRNA_Type,export_all_associations,run_from_scratch,force)
    
    null = EnsemblImport.getEnsemblAssociations(species,data_type,test_status); null=[]
    reformatExonFile(species,'exon',True); reformatExonFile(species,'junction',True)
    exportKnownJunctionComparisons(species)
    getExonAndJunctionSequences(species)

######### Below code deals with user read alignment as opposed to building the AltDatabase #########

class ExonInfo:
    def __init__(self,start,unique_id,annotation):
        self.start = start; self.unique_id = unique_id; self.annotation = annotation
    def ReadStart(self): return self.start
    def UniqueID(self): return self.unique_id
    def Annotation(self): return self.annotation
    def setExonRegionData(self,rd): self.rd = rd
    def ExonRegionData(self): return self.rd
    def setExonRegionID(self,region_id): self.region_id = region_id
    def ExonRegionID(self): return self.region_id
    def setAlignmentRegion(self,region_type): self.region_type = region_type
    def AlignmentRegion(self): return self.region_type
    def __repr__(self): return "ExonData values"
    
class JunctionData:
    def __init__(self,chr,strand,exon1_stop,exon2_start,junction_id,biotype):
        self.chr = chr; self.strand = strand; self._chr = chr
        self.exon1_stop = exon1_stop; self.exon2_start = exon2_start
        self.junction_id = junction_id; self.biotype = biotype
        #self.reads = reads; self.condition = condition
        self.left_exon = None; self.right_exon = None; self.jd = None; self.gene_id = None
        self.trans_splicing = None
        self.splice_events=''
        self.splice_junctions=''
        self.seq_length=''
        self.uid = None
    def Chr(self): return self.chr
    def Strand(self): return self.strand
    def Exon1Stop(self): return self.exon1_stop
    def Exon2Start(self): return self.exon2_start
    def setExon1Stop(self,exon1_stop): self.exon1_stop = exon1_stop
    def setExon2Start(self,exon2_start): self.exon2_start = exon2_start
    def setSeqLength(self,seq_length): self.seq_length = seq_length
    def SeqLength(self): return self.seq_length
    def BioType(self): return self.biotype
    def checkExonPosition(self,exon_pos):
        if exon_pos == self.Exon1Stop(): return 'left'
        else: return 'right'
    ### These are used to report novel exon boundaries
    def setExon1Start(self,exon1_start): self.exon1_start = exon1_start
    def setExon2Stop(self,exon2_stop): self.exon2_stop = exon2_stop
    def Exon1Start(self): return self.exon1_start
    def Exon2Stop(self): return self.exon2_stop
    def Reads(self): return self.reads
    def JunctionID(self): return self.junction_id
    def Condition(self): return self.condition
    def setExonAnnotations(self,jd):
        self.jd = jd
        self.splice_events = jd.AssociatedSplicingEvent()
        self.splice_junctions = jd.AssociatedSplicingJunctions()
        self.exon_region = jd.ExonRegionIDs()
        self.exonid = jd.ExonID()
        self.gene_id = jd.GeneID()
        self.uid = jd.GeneID()+':'+jd.ExonRegionIDs()
    def ExonAnnotations(self): return self.jd
    def setLeftExonAnnotations(self,ld): self.gene_id,self.left_exon = ld
    def LeftExonAnnotations(self): return self.left_exon
    def setRightExonAnnotations(self,rd): self.secondary_geneid,self.right_exon = rd
    def RightExonAnnotations(self): return self.right_exon
    def setGeneID(self,geneid): self.gene_id = geneid
    def GeneID(self): return self.gene_id
    def setSecondaryGeneID(self,secondary_geneid): self.secondary_geneid = secondary_geneid
    def SecondaryGeneID(self): return self.secondary_geneid
    def setTransSplicing(self): self.trans_splicing = 'yes'
    def TransSplicing(self): return self.trans_splicing
    def SpliceSitesFound(self):
        if self.jd != None: sites_found = 'both'
        elif self.left_exon != None and self.right_exon != None: sites_found = 'both'
        elif self.left_exon != None: sites_found = 'left'
        elif self.right_exon != None: sites_found = 'right'
        else: sites_found = None
        return sites_found
    def setConstitutive(self,constitutive): self.constitutive = constitutive
    def Constitutive(self): return self.constitutive
    def setAssociatedSplicingEvent(self,splice_events): self.splice_events = splice_events
    def AssociatedSplicingEvent(self): return self.splice_events
    def setAssociatedSplicingJunctions(self,splice_junctions): self.splice_junctions = splice_junctions
    def AssociatedSplicingJunctions(self): return self.splice_junctions
    def setExonID(self,exonid): self.exonid = exonid
    def ExonID(self): return self.exonid
    def setExonRegionID(self,exon_region): self.exon_region = exon_region
    def ExonRegionID(self): return self.exon_region
    def setUniqueID(self,uid): self.uid = uid
    def UniqueID(self): return self.uid
    def setLeftExonRegionData(self,li): self.li = li
    def LeftExonRegionData(self): return self.li
    def setRightExonRegionData(self,ri): self.ri = ri
    def RightExonRegionData(self): return self.ri
    def setNovel(self, side): self.side = side
    def Novel(self): return self.side    
    def __repr__(self): return "JunctionData values"

def checkBEDFileFormat(bed_dir,root_dir):
    """ This method checks to see if the BED files (junction or exon) have 'chr' proceeding the chr number.
    It also checks to see if some files have two underscores and one has none or if double underscores are missing from all."""
    dir_list = read_directory(bed_dir)
    x=0
    break_now = False
    chr_present = False
    condition_db={}
    for filename in dir_list:
        fn=filepath(bed_dir+filename)
        #if ('.bed' in fn or '.BED' in fn): delim = 'r'
        delim = 'rU'
        if '.tab' in string.lower(filename) or '.bed' in string.lower(filename):
            condition_db[filename]=[]
            for line in open(fn,delim).xreadlines(): ### changed rU to r to remove \r effectively, rather than read as end-lines
                if line[0] == '#': x=0 ### BioScope
                elif x == 0: x=1 ###skip the first line
                elif x < 10: ### Only check the first 10 lines
                    if 'chr' in line: ### Need to look at multiple input formats (chr could be in t[0] or t[1])
                        chr_present = True
                    x+=1
                else:
                    break_now = True
                    break
            if break_now == True:
                break
    
    ### Check to see if exon.bed and junction.bed file names are propper or faulty (which will result in downstream errors)
    double_underscores=[]
    no_doubles=[]
    for condition in condition_db:
        if '__' in condition:
            double_underscores.append(condition)
        else:
            no_doubles.append(condition)
    
    exon_beds=[]
    junctions_beds=[] 
    if len(double_underscores)>0 and len(no_doubles)>0:
        ### Hence, a problem is likely due to inconsistent naming
        print 'The input files appear to have inconsistent naming. If both exon and junction sample data are present, make sure they are named propperly.'
        print 'For example: cancer1__exon.bed, cancer1__junction.bed (double underscore required to match these samples up)!'
        print 'Exiting AltAnalyze'; forceError
    elif len(no_doubles)>0:
        for condition in no_doubles:
            condition = string.lower(condition)
            if 'exon' in condition:
                exon_beds.append(condition)
            if 'junction' in condition:
                junctions_beds.append(condition)
        if len(exon_beds)>0 and len(junctions_beds)>0:
            print 'The input files appear to have inconsistent naming. If both exon and junction sample data are present, make sure they are named propperly.'
            print 'For example: cancer1__exon.bed, cancer1__junction.bed (double underscore required to match these samples up)!'
            print 'Exiting AltAnalyze'; forceError
    return chr_present

def importBEDFile(bed_dir,root_dir,species,normalize_feature_exp,getReads=False,searchChr=None,getBiotype=None,testImport=False):
    dir_list = read_directory(bed_dir)
    begin_time = time.time()
    
    condition_count_db={}; neg_count=0; pos_count=0; junction_db={}; biotypes={}; algorithms={}; exon_len_db={}

    if testImport == 'yes': print "Reading user RNA-seq input data files"
    for filename in dir_list:
        count_db={}; x=0
        fn=filepath(bed_dir+filename)
        condition = export.findFilename(fn)
        if '__' in condition:
            ### Allow multiple junction files per sample to be combined (e.g. canonical and non-canonical junction alignments)
            condition=string.split(condition,'__')[0]+filename[-4:]
        if ('.bed' in fn or '.BED' in fn or '.tab' in fn or '.TAB' in fn) and '._' not in condition:
            if ('.bed' in fn or '.BED' in fn): delim = 'r'
            else: delim = 'rU'
            ### The below code removes .txt if still in the filename along with .tab or .bed
            if '.tab' in fn: condition = string.replace(condition,'.txt','.tab')
            elif '.bed' in fn: condition = string.replace(condition,'.txt','.bed')
            if '.TAB' in fn: condition = string.replace(condition,'.txt','.TAB')
            elif '.BED' in fn: condition = string.replace(condition,'.txt','.BED')
            
            if testImport == 'yes': print "Reading the bed file", [fn], condition
            ### If the BED was manually created on a Mac, will neeed 'rU' - test this
            for line in open(fn,delim).xreadlines(): break
            if len(line)>500: delim = 'rU'
            for line in open(fn,delim).xreadlines(): ### changed rU to r to remove \r effectively, rather than read as end-lines
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                if x==0 or '#' == data[0]:
                    format_description = data; x=1
                    algorithm = 'Unknown'
                    if 'TopHat' in format_description: algorithm = 'TopHat'
                    elif 'HMMSplicer' in format_description: algorithm = 'HMMSplicer'
                    elif 'SpliceMap junctions' in format_description: algorithm = 'SpliceMap'
                    elif t[0] == 'E1': algorithm = 'BioScope-junction'
                    elif '# filterOrphanedMates=' in data or 'alignmentFilteringMode=' in data or '#number_of_mapped_reads=' in data:
                        algorithm = 'BioScope-exon'
                    if testImport == 'yes': print condition, algorithm
                else:
                    try:
                        if searchChr == t[0] or ('BioScope' in algorithm and searchChr == t[1]): proceed = 'yes'
                        elif searchChr == 'chrMT':
                            if 'M' in t[0]: proceed = 'yes'
                            else: proceed = 'no'
                        else: proceed = 'no'
                    except IndexError:
                        print 'The input file:\n',filename
                        print 'is not formated as expected (format='+algorithm+').'
                        print 'search chromosome:',searchChr
                        print t; force_bad_exit
                    if proceed == 'yes':
                        proceed = 'no'
                        if '.tab' in fn or '.TAB' in fn:
                            ### Applies to non-BED format Junction and Exon inputs (BioScope)
                            if 'BioScope' in algorithm:
                                if algorithm == 'BioScope-exon': ### Not BED format
                                    chr,source,data_type,start,end,reads,strand,null,gene_info=t[:9]
                                    if data_type == 'exon': ### Can also be CDS
                                        gene_info,test,rpkm_info,null = string.split(gene_info,';')
                                        symbol = string.split(gene_info,' ')[-1]
                                        #refseq = string.split(transcript_info,' ')[-1]
                                        rpkm = string.split(rpkm_info,' ')[-1]
                                        #if normalize_feature_exp == 'RPKM': reads = rpkm ### The RPKM should be adjusted +1 counts, so don't use this
                                        biotype = 'exon'; biotypes[biotype]=[]
                                        exon1_stop,exon2_start = int(start),int(end); junction_id=''
                                        ### Adjust exon positions - not ideal but necessary. Needed as a result of exon regions overlapping by 1nt (due to build process)
                                        exon1_stop+=1; exon2_start-=1
                                        if float(reads)>0: proceed = 'yes'
                                        seq_length = abs(exon1_stop-exon2_start)
                                if algorithm == 'BioScope-junction':
                                    chr = t[1]; strand = t[2]; exon1_stop = int(t[4]); exon2_start = int(t[8]); count_paired = t[17]; count_single = t[19]; score=t[21]
                                    try: exon1_start = int(t[3]); exon2_stop = int(t[9])
                                    except Exception: null=[] ### If missing, these are not assigned
                                    reads = str(int(float(count_paired))+int(float(count_single))) ### Users will either have paired or single read (this uses either)
                                    biotype = 'junction'; biotypes[biotype]=[]; junction_id=''
                                    if float(reads)>0: proceed = 'yes'
                                    seq_length = abs(float(exon1_stop-exon2_start))
                        else:
                            try:
                                ### Applies to BED format Junction input
                                chr, exon1_start, exon2_stop, junction_id, reads, strand, null, null, null, null, lengths, null = t
                                exon1_len,exon2_len=string.split(lengths,',')[:2]; exon1_len = int(exon1_len); exon2_len = int(exon2_len)
                                exon1_start = int(exon1_start); exon2_stop = int(exon2_stop)
                                biotype = 'junction'; biotypes[biotype]=[]
                                if strand == '-':
                                    exon1_stop = exon1_start+exon1_len; exon2_start=exon2_stop-exon2_len+1
                                    ### Exons have the opposite order
                                    a = exon1_start,exon1_stop; b = exon2_start,exon2_stop
                                    exon1_stop,exon1_start = b; exon2_stop,exon2_start = a
                                else:
                                    exon1_stop = exon1_start+exon1_len; exon2_start=exon2_stop-exon2_len+1
                                proceed = 'yes'
                                if algorithm == 'HMMSplicer':
                                    if '|junc=' in junction_id: reads = string.split(junction_id,'|junc=')[-1]
                                    else: proceed = 'no'
                                if algorithm == 'SpliceMap':
                                    if ')' in junction_id and len(junction_id)>1: reads = string.split(junction_id,')')[0][1:]
                                    else: proceed = 'no'
                                seq_length = abs(float(exon1_stop-exon2_start)) ### Junction distance
                            except Exception:
                                ### Applies to BED format exon input (BEDTools export)
                                # bamToBed -i accepted_hits.bam -split| coverageBed -a stdin -b /home/nsalomonis/databases/Mm_Ensembl_exons.bed > day0_8B__exons.bed
                                try: chr, start, end, exon_id, null, strand, reads, bp_coverage, bp_total, percent_coverage = t
                                except Exception:
                                    print 'The file',fn,'does not appear to be propperly formatted as input.'
                                    print t; force_exception
                                algorithm = 'TopHat-exon'; biotype = 'exon'; biotypes[biotype]=[]
                                exon1_stop,exon2_start = int(start),int(end); junction_id=exon_id; seq_length = float(bp_total)
                                if seq_length == 0:
                                    seq_length = abs(float(exon1_stop-exon2_start))
                                ### Adjust exon positions - not ideal but necessary. Needed as a result of exon regions overlapping by 1nt (due to build process)
                                exon1_stop+=1; exon2_start-=1
                                if float(reads)>0:
                                    proceed = 'yes'
                                else: proceed = 'no'

                        if proceed == 'yes':
                            if 'chr' not in chr:
                                chr = 'chr'+chr ### Add the chromosome prefix
                            if chr == 'chrM': chr = 'chrMT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
                            if chr == 'M': chr = 'MT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
                            if strand == '+': pos_count+=1
                            else: neg_count+=1
                            if getReads and seq_length>0:
                                if getBiotype == biotype:
                                    count_db[chr,exon1_stop,exon2_start] = reads
                                    try: exon_len_db[chr,exon1_stop,exon2_start] = seq_length
                                    except Exception: exon_len_db[chr,exon1_stop,exon2_start] = []
                            elif seq_length>0:
                                if (chr,exon1_stop,exon2_start) not in junction_db:
                                    ji = JunctionData(chr,strand,exon1_stop,exon2_start,junction_id,biotype)
                                    junction_db[chr,exon1_stop,exon2_start] = ji
                                    try: ji.setSeqLength(seq_length) ### If RPKM imported or calculated
                                    except Exception: null=[]
                                    try: ji.setExon1Start(exon1_start);ji.setExon2Stop(exon2_stop)
                                    except Exception: null=[]
                                    key = chr,exon1_stop,exon2_start
                algorithms[algorithm]=[]
            if getReads:
                if condition in condition_count_db:
                    ### combine the data from the different files for the same sample junction alignments
                    count_db1 = condition_count_db[condition]
                    for key in count_db:
                        if key not in count_db1: count_db1[key] = count_db[key]
                        else:
                            combined_counts = int(count_db1[key])+int(count_db[key])
                            count_db1[key] = str(combined_counts)
                    condition_count_db[condition]=count_db1
                else:
                    try: condition_count_db[condition] = count_db
                    except Exception: null=[] ### Occurs for other text files in the directory that are not used for the analysis
                
    end_time = time.time()
    if testImport == 'yes': print 'Read coordinates imported in',int(end_time-begin_time),'seconds'
    if getReads:
        #print len(exon_len_db), getBiotype, 'read counts present for',algorithm
        return condition_count_db,exon_len_db,biotypes,algorithms
    else:
        if testImport == 'yes':
            if 'exon' not in biotypes and 'BioScope' not in algorithm:
                print len(junction_db),'junctions present in',algorithm,'format BED files.' # ('+str(pos_count),str(neg_count)+' by strand).'
            elif 'exon' in biotypes and 'BioScope' not in algorithm:
                print len(junction_db),'sequence identifiers present in input files.' 
            else: print len(junction_db),'sequence identifiers present in BioScope input files.'
        return junction_db,biotypes,algorithms

def importExonCoordinates(probeCoordinateFile,search_chr,getBiotype):
    probe_coordinate_db={}
    junction_db={}
    biotypes={}
    x=0
    fn=filepath(probeCoordinateFile)
    for line in open(fn,'rU').xreadlines(): ### changed rU to r to remove \r effectively, rather than read as end-lines
        data = cleanUpLine(line)
        if x==0: x=1
        else:
            t = string.split(data,'\t')
            probe_id = t[0]; probeset_id=t[1]; chr=t[2]; strand=t[3]; start=t[4]; end=t[5]
            exon1_stop,exon2_start = int(start),int(end)
            seq_length = abs(float(exon1_stop-exon2_start))
            if 'chr' not in chr:
                chr = 'chr'+chr ### Add the chromosome prefix
            if chr == 'chrM': chr = 'chrMT' ### MT is the Ensembl convention whereas M is the Affymetrix and UCSC convention
            if search_chr == chr or search_chr == None:
                try: biotype = t[6]
                except Exception:
                    if seq_length>25:biotype = 'junction'
                    else: biotype = 'exon'
                if strand == '-':
                    exon1_stop,exon2_start = exon2_start, exon1_stop ### this is their actual 5' -> 3' orientation
                if biotype == 'junction':
                    exon1_start,exon2_stop = exon1_stop,exon2_start
                else:
                    exon1_stop+=1; exon2_start-=1
                biotypes[biotype]=[]
                if getBiotype == biotype or getBiotype == None:
                    ji = JunctionData(chr,strand,exon1_stop,exon2_start,probe_id,biotype)
                    junction_db[chr,exon1_stop,exon2_start] = ji
                    try: ji.setSeqLength(seq_length) ### If RPKM imported or calculated
                    except Exception: null=[]
                    try: ji.setExon1Start(exon1_start);ji.setExon2Stop(exon2_stop)
                    except Exception: null=[]
                    probe_coordinate_db[probe_id] = chr,exon1_stop,exon2_start ### Import the expression data for the correct chromosomes with these IDs

    return probe_coordinate_db, junction_db, biotypes
       
def importExpressionMatrix(exp_dir,root_dir,species,fl,getReads,search_chr=None,getBiotype=None):
    """ Non-RNA-Seq expression data (typically Affymetrix microarray) import and mapping to an external probe-coordinate database """
    begin_time = time.time()
            
    condition_count_db={}; neg_count=0; pos_count=0; algorithms={}; exon_len_db={}
    
    probe_coordinate_db, junction_db, biotypes = importExonCoordinates(fl.ExonMapFile(),search_chr,getBiotype)
    
    x=0
    fn=filepath(exp_dir)[:-1]
    condition = export.findFilename(fn)
    ### If the BED was manually created on a Mac, will neeed 'rU' - test this
    for line in open(fn,'rU').xreadlines(): ### changed rU to r to remove \r effectively, rather than read as end-lines
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if '#' == data[0]: None
        elif x==0:
            if 'block' in t:
                start_index = 7
            else:
                start_index = 1
            headers = t[start_index:]
            x=1
        else:
            proceed = 'yes' ### restrict by chromosome with minimum line parsing (unless we want counts instead)
            probe_id=t[0]
            if probe_id in probe_coordinate_db:
                key = probe_coordinate_db[probe_id]
                if getReads == 'no':
                    pass
                else:
                    expression_data = t[start_index:]
                    i=0
                    for sample in headers:
                        if sample in condition_count_db:
                            count_db = condition_count_db[sample]
                            count_db[key] = expression_data[i]
                            exon_len_db[key]=[]
                        else:
                            count_db={}
                            count_db[key] = expression_data[i]
                            condition_count_db[sample] = count_db
                            exon_len_db[key]=[]
                        i+=1

    algorithms['ProbeData']=[]
    end_time = time.time()
    if testImport == 'yes': print 'Probe data imported in',int(end_time-begin_time),'seconds'
    if getReads == 'yes':
        return condition_count_db,exon_len_db,biotypes,algorithms
    else:
        return junction_db,biotypes,algorithms

def adjustCounts(condition_count_db,exon_len_db):
    for key in exon_len_db:
        try:
            null=exon_len_db[key]
            for condition in condition_count_db:
                count_db = condition_count_db[condition]
                try: read_count = float(count_db[key])+1 ###This adjustment allows us to obtain more realist folds where 0 is compared and use log2
                except KeyError: read_count = 1 ###Was zero, but needs to be one for more realistic log2 fold calculations
                count_db[key] = str(read_count) ### Replace original counts with adjusted counts
        except Exception: null=[]
    return condition_count_db

def calculateRPKM(condition_count_db,exon_len_db,biotype_to_examine):
    """Determines the total number of reads in a sample and then calculates RPMK relative to a pre-determined junction length (60).
    60 was choosen, based on Illumina single-end read lengths of 35 (5 nt allowed overhand on either side of the junction)"""
    ### Get the total number of mapped reads
    mapped_reads={}
    for condition in condition_count_db:
        mapped_reads[condition]=0
        count_db = condition_count_db[condition]
        for key in count_db:
            read_count = count_db[key]
            mapped_reads[condition]+=float(read_count)
            
    ### Use the average_total_reads when no counts reported such that 0 counts are comparable            
    average_total_reads = 0
    for i in mapped_reads:
        average_total_reads+=mapped_reads[i]
        if testImport == 'yes':
            print 'condition:',i,'total reads:',mapped_reads[i]
    average_total_reads = average_total_reads/len(condition_count_db)
    if testImport == 'yes':
        print 'average_total_reads:',average_total_reads
    k=0
    c=math.pow(10.0,9.0)
    for key in exon_len_db:
        try:
            for condition in condition_count_db:
                total_mapped_reads = mapped_reads[condition]
                try: read_count = float(condition_count_db[condition][key])+1 ###This adjustment allows us to obtain more realist folds where 0 is compared and use log2
                except KeyError: read_count = 1 ###Was zero, but needs to be one for more realistic log2 fold calculations
                if biotype_to_examine == 'junction': region_length = 60.0
                else:
                    try: region_length = exon_len_db[key]
                    except Exception: continue ### This should only occur during testing (when restricting to one or few chromosomes)
                if read_count == 1: ###This adjustment allows us to obtain more realist folds where 0 is compared and use log2
                    rpkm = c*(float(read_count)/(float(average_total_reads)*region_length)) 
                try:
                    if region_length == 0:
                        region_length = abs(int(key[2]-key[1]))
                    rpkm = c*(read_count/(float(total_mapped_reads)*region_length))
                except Exception:
                    print condition, key
                    print 'Error Encountered... Exon or Junction of zero length encoutered... RPKM failed... Exiting AltAnalyze.'
                    print 'This error may be due to inconsistent file naming. If both exon and junction sample data is present, make sure they are named propperly.'
                    print 'For example: cancer1__exon.bed, cancer1__junction.bed (double underscore required to match these samples up)!'
                    print [read_count,total_mapped_reads,region_length];k=1; forceError
                condition_count_db[condition][key] = str(rpkm) ### Replace original counts with RPMK
        except Exception:
            if k == 1: kill
            null=[]
    return condition_count_db

def calculateGeneRPKM(gene_count_db):
    """Determines the total number of reads in a sample and then calculates RPMK relative to a pre-determined junction length (60).
    60 was choosen, based on Illumina single-end read lengths of 35 (5 nt allowed overhand on either side of the junction)"""
    ### Get the total number of mapped reads (relative to all gene aligned rather than genome aligned exon reads)
    mapped_reads={}
    for gene in gene_count_db:
        index=0
        for (read_count,total_len) in gene_count_db[gene]:
            try: mapped_reads[index]+=float(read_count)
            except Exception: mapped_reads[index]=float(read_count)
            index+=1
            
    ### Use the average_total_reads when no counts reported such that 0 counts are comparable            
    average_total_reads = 0
    for i in mapped_reads: average_total_reads+=mapped_reads[i]
    average_total_reads = average_total_reads/(index+1)
    
    c=math.pow(10.0,9.0)
    for gene in gene_count_db:
        index=0; rpkms = []
        for (read_count,region_length) in gene_count_db[gene]:
            total_mapped_reads = mapped_reads[index]
            #print [read_count],[region_length],[total_mapped_reads]
            #if gene == 'ENSMUSG00000028186': print [read_count, index, total_mapped_reads,average_total_reads,region_length]
            if read_count == 0: read_count=1; rpkm = c*(float(read_count)/(float(average_total_reads)*region_length)) ###This adjustment allows us to obtain more realist folds where 0 is compared and use log2
            else:
                try: rpkm = c*(float(read_count+1)/(float(total_mapped_reads)*region_length)) ### read count is incremented +1 (see next line)
                except Exception: read_count=1; rpkm = c*(float(read_count)/(float(average_total_reads)*region_length)) ###This adjustment allows us to obtain more realist folds where 0 is compared and use log2
            #if gene == 'ENSMUSG00000028186': print rpkm,read_count,index,total_mapped_reads,average_total_reads,region_length
            rpkms.append(rpkm)
            index+=1
        gene_count_db[gene] = rpkms ### Replace original counts with RPMK
    return gene_count_db

def calculateGeneLevelStatistics(steady_state_export,expressed_gene_exon_db,normalize_feature_exp,array_names,fl):
    global UserOptions; UserOptions = fl
    exp_file = string.replace(steady_state_export,'-steady-state','')
    if normalize_feature_exp == 'RPKM':
        exp_dbase, all_exp_features, array_count = importRawCountData(exp_file,expressed_gene_exon_db)
        steady_state_db = obtainGeneCounts(expressed_gene_exon_db,exp_dbase,array_count,normalize_feature_exp); exp_dbase=[]
        exportGeneCounts(steady_state_export,array_names,steady_state_db)
        steady_state_db = calculateGeneRPKM(steady_state_db)
    else:
        exp_dbase, all_exp_features, array_count = importNormalizedCountData(exp_file,expressed_gene_exon_db)
        steady_state_db = obtainGeneCounts(expressed_gene_exon_db,exp_dbase,array_count,normalize_feature_exp); exp_dbase=[]
        exportGeneCounts(steady_state_export,array_names,steady_state_db)
    return steady_state_db, all_exp_features
    
def exportGeneCounts(steady_state_export,headers,gene_count_db):
    ### In addition to RPKM gene-level data, export gene level counts and lengths (should be able to calculate gene RPKMs from this file)
    export_path = string.replace(steady_state_export,'exp.','counts.')
    export_data = export.ExportFile(export_path)
    
    title = string.join(['Ensembl']+headers,'\t')+'\n'
    export_data.write(title)
        
    for gene in gene_count_db:
        sample_counts=[]
        for count_data in gene_count_db[gene]:
            try: read_count,region_length = count_data
            except Exception: read_count = count_data
            sample_counts.append(str(read_count))
        sample_counts = string.join([gene]+sample_counts,'\t')+'\n'
        export_data.write(sample_counts)
    export_data.close()

def importGeneCounts(filename,import_type):
    ### Import non-normalized original counts and return the max value            
    counts_filename = string.replace(filename,'exp.','counts.')
    status = verifyFile(counts_filename)
    if status == 'not found': ### Occurs for non-normalized counts
        counts_filename = filename
    fn=filepath(counts_filename); x=0; count_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: array_names = t[1:]; x=1
        else:
            gene = t[0]
            if import_type == 'max':
                count_db[gene] = str(max(map(float,t[1:])))
            else:
                count_db[gene] = map(float,t[1:])
    return count_db,array_names

def calculateGeneRPKM(gene_count_db):
    """Determines the total number of reads in a sample and then calculates RPMK relative to a pre-determined junction length (60).
    60 was choosen, based on Illumina single-end read lengths of 35 (5 nt allowed overhand on either side of the junction)"""
    ### Get the total number of mapped reads (relative to all gene aligned rather than genome aligned exon reads)
    mapped_reads={}
    for gene in gene_count_db:
        index=0
        for (read_count,total_len) in gene_count_db[gene]:
            try: mapped_reads[index]+=float(read_count)
            except Exception: mapped_reads[index]=float(read_count)
            index+=1
            
    ### Use the average_total_reads when no counts reported such that 0 counts are comparable            
    average_total_reads = 0
    for i in mapped_reads: average_total_reads+=mapped_reads[i]
    average_total_reads = average_total_reads/(index+1) ###
    
    c=math.pow(10.0,9.0)
    for gene in gene_count_db:
        index=0; rpkms = []
        for (read_count,region_length) in gene_count_db[gene]:
            total_mapped_reads = mapped_reads[index]
            #print [read_count],[region_length],[total_mapped_reads]
            #if gene == 'ENSMUSG00000028186': print [read_count, index, total_mapped_reads,average_total_reads,region_length]
            if read_count == 0: read_count=1; rpkm = c*(float(read_count)/(float(average_total_reads)*region_length)) ###This adjustment allows us to obtain more realist folds where 0 is compared and use log2
            else:
                try: rpkm = c*(float(read_count+1)/(float(total_mapped_reads)*region_length)) ### read count is incremented +1 (see next line)
                except Exception: read_count=1; rpkm = c*(float(read_count)/(float(average_total_reads)*region_length)) ###This adjustment allows us to obtain more realist folds where 0 is compared and use log2
            #if gene == 'ENSMUSG00000028186': print rpkm,read_count,index,total_mapped_reads,average_total_reads,region_length
            rpkms.append(rpkm)
            index+=1
        gene_count_db[gene] = rpkms ### Replace original counts with RPMK
    return gene_count_db

def deleteOldAnnotations(species,root_dir,dataset_name):
    db_dir = root_dir+'AltDatabase/'+species
    try:
        status = export.deleteFolder(db_dir)
        if status == 'success':
            print "...Previous experiment database deleted"
    except Exception: null=[]
    
    count_dir = root_dir+'ExpressionInput/Counts'
    try: status = export.deleteFolder(count_dir)
    except Exception: pass
    
    if 'exp.' not in dataset_name: dataset_name = 'exp.'+dataset_name
    if '.txt' not in dataset_name: dataset_name+='.txt'
    export_path = root_dir+'ExpressionInput/'+dataset_name
    try: os.remove(filepath(export_path))
    except Exception: null=[]
    try: os.remove(filepath(string.replace(export_path,'exp.','counts.')))
    except Exception: null=[]
    try: os.remove(filepath(string.replace(export_path,'exp.','novel.')))
    except Exception: null=[]

from copy_reg import pickle
from types import MethodType

def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)

def call_it(instance, name, args=(), kwargs=None):
    "indirect caller for instance methods and multiprocessing"
    if kwargs is None:
        kwargs = {}
    return getattr(instance, name)(*args, **kwargs)

def alignExonsAndJunctionsToEnsembl(species,exp_file_location_db,dataset_name,Multi=None):
    fl = exp_file_location_db[dataset_name]
    
    try: multiThreading = fl.multiThreading()
    except Exception: multiThreading = True
    print 'multiThreading:',multiThreading
    
    normalize_feature_exp = fl.FeatureNormalization()
    
    testImport='no'
    rnaseq_begin_time = time.time()
    p = AlignExonsAndJunctionsToEnsembl(species,exp_file_location_db,dataset_name,testImport)
    chromosomes = p.getChromosomes()
    
    ### The following files need to be produced from chromosome specific sets later
    countsFile = p.countsFile()
    exonFile = p.exonFile()
    junctionFile = p.junctionFile()
    junctionCompFile = p.junctionCompFile()
    
    #chromosomes = ['chrY','chr1','chr2']
    #p('chrY'); p('chr1'); p('chr2')
    #chromosomes = ['chr8','chr17']
    multiprocessing_pipe = True
    
    if 'exp.' not in dataset_name:
        dataset_name = 'exp.'+dataset_name
        if '.txt' not in dataset_name:
            dataset_name+='.txt'
                
    try:
        mlp=Multi
        pool_size = mlp.cpu_count()
        if multiprocessing_pipe and multiThreading:
            ### This is like pool, but less efficient (needed to get print outs)
            s = pool_size; b=0
            chr_blocks=[]
            while s<len(chromosomes):
                chr_blocks.append(chromosomes[b:s])
                b+=pool_size; s+=pool_size
            chr_blocks.append(chromosomes[b:s])
    
            queue = mlp.Queue()
            results=[]
            #parent_conn, child_conn=multiprocessing.Pipe()
            for chromosomes in chr_blocks:
                procs=list()
                #print 'Block size:',len(chromosomes)
                for search_chr in chromosomes:
                    proc = mlp.Process(target=p, args=(queue,search_chr)) ### passing sys.stdout unfortunately doesn't work to pass the Tk string 
                    procs.append(proc)
                    proc.start()

                for _ in procs:
                    val = queue.get()
                    if p.AnalysisMode() == 'GUI': print '*',
                    results.append(val)
                    
                for proc in procs:
                    proc.join()
        
        elif multiThreading:
            pool = mlp.Pool(processes=pool_size)
            chr_vars=[]
            for search_chr in chromosomes:
                chr_vars.append(([],search_chr)) ### As an alternative for the pipe version above, pass an empty list rather than queue
            
            results = pool.map(p, chr_vars) ### worker jobs initiated in tandem
            try:pool.close(); pool.join(); pool = None
            except Exception: pass
        else:
            forceThreadingError
        print 'Read exon and junction mapping complete'
        
    except Exception,e:
        #print e
        print 'Proceeding with single-processor version align...'
        try: proc.close; proc.join; proc = None
        except Exception: pass
        try: pool.close(); pool.join(); pool = None
        except Exception: pass
        results=[] ### For single-thread compatible versions of Python
        for search_chr in chromosomes:
            result = p([],search_chr)
            results.append(result)
        
    results_organized=[]
    for result_set in results:
        if len(result_set[0])>0: ### Sometimes chromsomes are missing
            biotypes = result_set[0]
        results_organized.append(list(result_set[1:]))
    
    pooled_results = [sum(value) for value in zip(*results_organized)] # combine these counts
    pooled_results = [biotypes]+pooled_results
    p.setCountsOverview(pooled_results) # store as retreivable objects

    catFiles(countsFile,'Counts')
    catFiles(junctionFile,'junctions')
    catFiles(exonFile,'exons')
    catFiles(junctionCompFile,'comps')
    
    if normalize_feature_exp == 'RPKM':
        fastRPKMCalculate(countsFile)
    
    rnaseq_end_time = time.time()
    print '...RNA-seq import completed in',int(rnaseq_end_time-rnaseq_begin_time),'seconds\n'
    biotypes = p.outputResults()
    return biotypes

def catFiles(outFileDir,folder):
    """ Concatenate all the chromosomal files but retain only the first header """
    root_dir = export.findParentDir(outFileDir)+folder+'/'
    dir_list = read_directory(root_dir)
    firstFile=True
    with open(filepath(outFileDir), 'w') as outfile:
        for fname in dir_list:
            chr_file = root_dir+fname
            header=True
            with open(filepath(chr_file)) as infile:
                for line in infile:
                    if header:
                        header=False
                        if firstFile:
                            outfile.write(line)
                            firstFile=False
                    else: outfile.write(line)
    export.deleteFolder(root_dir)
                
def error(msg, *args):
    return multiprocessing.get_logger().error(msg, *args)

class AlignExonsAndJunctionsToEnsembl:
    def setCountsOverview(self, overview):
        self.biotypes_store, self.known_count, self.novel_junction_count, self.trans_splicing_reads, self.junctions_without_exon_gene_alignments, self.exons_without_gene_alignment_count = overview
    
    def getChromosomes(self):
        chr_list=list()
        for c in self.chromosomes:
            ### Sort chromosome by int number
            ci=string.replace(c,'chr','')
            try: ci = int(ci)
            except Exception: pass
            chr_list.append((ci,c))
        chr_list.sort()
        chr_list2=list()
        for (i,c) in chr_list: chr_list2.append(c) ### sorted
        return chr_list2
    
    def countsFile(self):
        return string.replace(self.expfile,'exp.','counts.')
    
    def junctionFile(self):
        junction_file = self.root_dir+'AltDatabase/'+self.species+'/RNASeq/'+self.species + '_Ensembl_junctions.txt'
        return junction_file

    def exonFile(self):
        exon_file = self.root_dir+'AltDatabase/'+self.species+'/RNASeq/'+self.species + '_Ensembl_exons.txt'
        return exon_file

    def junctionCompFile(self):
        junction_comp_file = self.root_dir+'AltDatabase/'+self.species+'/RNASeq/'+self.species + '_junction_comps_updated.txt'
        return junction_comp_file
    
    def AnalysisMode(self): return self.analysisMode
    
    def __init__(self,species,exp_file_location_db,dataset_name,testImport):
        self.species = species; self.dataset_name = dataset_name
        self.testImport = testImport
        
        fl = exp_file_location_db[dataset_name]
        bed_dir=fl.BEDFileDir()
        root_dir=fl.RootDir()
        #self.stdout = fl.STDOUT()
        try: platformType = fl.PlatformType()
        except Exception: platformType = 'RNASeq'
        try: analysisMode = fl.AnalysisMode()
        except Exception: analysisMode = 'GUI'
        
        ### This occurs when run using the BAMtoBED pipeline in the GUI
        if 'exp.' not in dataset_name:
            dataset_name = 'exp.'+dataset_name
            if '.txt' not in dataset_name:
                dataset_name+='.txt'
            self.dataset_name = dataset_name
            
        ### Import experimentally identified junction splice-sites
        normalize_feature_exp = fl.FeatureNormalization()

        if platformType == 'RNASeq':
            chr_status = checkBEDFileFormat(bed_dir,root_dir) ### If false, need to remove 'chr' from the search_chr
        else:
            chr_status = True

        #self.fl = fl # Can not pass this object in pool or it breaks
        self.platformType = platformType
        self.analysisMode = analysisMode
        self.root_dir = root_dir
        self.normalize_feature_exp = normalize_feature_exp
        self.bed_dir = bed_dir
        self.chr_status = chr_status
        self.exonBedBuildStatus = fl.ExonBedBuildStatus()
        self.expfile = root_dir+'ExpressionInput/'+dataset_name
        
        if testImport == 'yes':
            print 'Chromosome annotation detected =',chr_status
        if self.exonBedBuildStatus == 'yes':
            reformatExonFile(species,'exon',chr_status) ### exports BED format exons for exon expression extraction
        """
        Strategies to reduce memory in RNASeq:
        1)	(done)Delete old AltDatabase-local version if it exists before starting
        2)	(done)Check to see if a file exists before writing it and if so append rather than create
        3)	(done)Get counts last and normalize last in for exons and junctions separately.
        4)	(done)Delete objects explicitly before importing any new data (define a new function that just does this).
        5)	(done)Get all chromosomes first then parse exon and junction coordinate data on a per known chromosome basis.
        6)	(done)Prior to deleting all junction/exon object info for each chromsome, save the coordinate(key)-to-annotation information for the read count export file."""
        
        ### Delete any existing annotation databases that currently exist (redundant with below)
        deleteOldAnnotations(species,root_dir,dataset_name)
        
        ###Define variables to report once reads for all chromosomes have been aligned
        #global self.known_count; global self.novel_junction_count; global self.one_found; global self.not_found; global self.both_found; global self.trans_splicing_reads
        #global self.junctions_without_exon_gene_alignments; global self.exons_without_gene_alignment_count; global self.junction_simple_db; global self.chr_strand_gene_dbs

        self.known_count=0; self.novel_junction_count=0; self.one_found=0; self.not_found=0; self.both_found=0; self.trans_splicing_reads=0
        self.junctions_without_exon_gene_alignments=0; self.exons_without_gene_alignment_count=0; self.junction_simple_db={}
        
        ###Begin Chromosome specific read to exon alignments
        self.chr_strand_gene_dbs,self.location_gene_db,chromosomes,self.gene_location_db = getChromosomeStrandCoordinates(species,testImport)
        self.chromosomes = chromosomes

        print "Processing exon/junction coordinates sequentially by chromosome"
        print "Note: this step is time intensive (can be hours) and no print statements may post for a while"
        
    def outputResults(self):
        exportDatasetLinkedGenes(self.species,self.gene_location_db,self.root_dir) ### Include an entry for gene IDs to include constitutive expression for RPKM normalized data     
        chr_gene_locations=[]; self.location_gene_db=[]; self.chr_strand_gene_dbs=[] 
        
        #print 'user coordinates imported/processed'
        #print 'Importing read counts from coordinate data...'

        biotypes = self.biotypes_store
        ### Output summary statistics
        if self.normalize_feature_exp != 'none':
            print self.normalize_feature_exp, 'normalization complete'
        if 'junction' in biotypes:
            print 'Imported Junction Statistics:'
            print '    ',self.known_count, 'junctions found in Ensembl/UCSC and',self.novel_junction_count,'are novel'
            print '    ',self.trans_splicing_reads,'trans-splicing junctions found (two aligning Ensembl genes)'
            print '    ',self.junctions_without_exon_gene_alignments, 'junctions where neither splice-site aligned to a gene'
            if (float(self.known_count)*10)<float(self.novel_junction_count):
                print '\nWARNING!!!!! Few junctions aligned to known exons. Ensure that the AltAnalyze Ensembl database\nversion matches the genome build aligned to!\n'
        if 'exon' in biotypes:
            print 'Imported Exon Statistics:'
            print '    ',self.exons_without_gene_alignment_count, 'exons where neither aligned to a gene'
        print 'User databases and read counts written to:', self.root_dir[:-1]+'ExpressionInput'
        
        
        ### END CHROMOSOME SPECIFIC ANALYSES
        if self.exonBedBuildStatus == 'yes':
            bedfile = moveBAMtoBEDFile(self.species,self.dataset_name,self.root_dir)
            print 'Exon BED file updated with novel exon predictions from junction file'
            return bedfile; sys.exit()
                
        clearObjectsFromMemory(self.junction_simple_db); self.junction_simple_db=[]
        return biotypes

    def test(self, search_chr):
        print search_chr
        
    def __call__(self, queue, search_chr):
        try:
            #sys.stdout = self.stdout
            platformType = self.platformType
            testImport = self.testImport
            species = self.species
            dataset_name = self.dataset_name
            
            platformType = self.platformType
            analysisMode = self.analysisMode
            root_dir = self.root_dir
            normalize_feature_exp = self.normalize_feature_exp
            bed_dir = self.bed_dir
            chr_status = self.chr_status
            
            junction_annotations={}
            if chr_status == False:
                searchchr = string.replace(search_chr,'chr','')
            else:
                searchchr = search_chr
            if platformType == 'RNASeq':
                junction_db,biotypes,algorithms = importBEDFile(bed_dir,root_dir,species,normalize_feature_exp,searchChr=searchchr,testImport=testImport)
            else:
                normalize_feature_exp = 'quantile'
                junction_db,biotypes,algorithms = importExpressionMatrix(bed_dir,root_dir,species,fl,'no',search_chr=searchchr)
            
            self.biotypes_store = biotypes
            if len(junction_db)>0:
                    
                ### Determine which kind of data is being imported, junctions, exons or both
                unmapped_exon_db={}
                    
                if 'junction' in biotypes:
                    ### Get all known junction splice-sites
                    ens_junction_coord_db = importExonAnnotations(species,'junction_coordinates',search_chr)
                    if testImport == 'yes':
                        print len(ens_junction_coord_db),'Ensembl/UCSC junctions imported'
                ### Identify known junctions sites found in the experimental dataset (perfect match)
                novel_junction_db={}; novel_exon_db={}
                for key in junction_db:
                    ji=junction_db[key]
                    if ji.BioType()=='junction':
                        if key in ens_junction_coord_db:
                            jd=ens_junction_coord_db[key]
                            ji.setExonAnnotations(jd)
                            self.known_count+=1
                        else:
                            novel_junction_db[key]=junction_db[key]; self.novel_junction_count+=1
                            #if 75953254 in key: print key; sys.exit()
                    else:
                        unmapped_exon_db[key]=junction_db[key]
                ens_exon_db = importExonAnnotations(species,'exon',search_chr)
                
                if 'junction' in biotypes:
                    if testImport == 'yes':
                        print self.known_count,  'junctions found in Ensembl/UCSC and',len(novel_junction_db),'are novel.'
                    
                    ### Separate each junction into a 5' and 3' splice site (exon1_coord_db and exon2_coord_db)
                    exon1_coord_db={}; exon2_coord_db={}    
                    for (chr,exon1_stop,exon2_start) in ens_junction_coord_db:
                        jd = ens_junction_coord_db[(chr,exon1_stop,exon2_start)]
                        exon1_coord_db[chr,exon1_stop] = jd.GeneID(),string.split(jd.ExonRegionIDs(),'-')[0]
                        exon2_coord_db[chr,exon2_start] = jd.GeneID(),string.split(jd.ExonRegionIDs(),'-')[1]
                    clearObjectsFromMemory(ens_junction_coord_db); ens_junction_coord_db=[] ### Clear object from memory
                    
                    ### Get and re-format individual exon info
                    exon_region_db={}
                    #if 'exon' not in biotypes:
                    for gene in ens_exon_db:
                        for rd in ens_exon_db[gene]:
                            exon_region_db[gene,rd.ExonRegionIDs()]=rd
        
                    ### Add the exon annotations from the known junctions to the exons to export dictionary
                    exons_to_export={}
                    for key in junction_db:
                        ji=junction_db[key]
                        if ji.ExonAnnotations() != None:
                            jd = ji.ExonAnnotations()
                            exon1, exon2 = string.split(jd.ExonRegionIDs(),'-')
                            key1 = jd.GeneID(),exon1; key2 = jd.GeneID(),exon2
                            exons_to_export[key1] = exon_region_db[key1]
                            exons_to_export[key2] = exon_region_db[key2]
                            
                    ### For novel experimental junctions, identify those with at least one matching known 5' or 3' site
                    exons_not_identified = {}; novel_exon_coordinates={}
                    for (chr,exon1_stop,exon2_start) in novel_junction_db:
                        ji = novel_junction_db[(chr,exon1_stop,exon2_start)]
                        coord = [exon1_stop,exon2_start]; coord.sort()
                        if (chr,exon1_stop) in exon1_coord_db and (chr,exon2_start) in exon2_coord_db:
                            ### Assign exon annotations to junctions where both splice-sites are known in Ensembl/UCSC
                            ### Store the exon objects, genes and regions (le is a tuple of gene and exon region ID)
                            ### Do this later for the below un-assigned exons
                            le=exon1_coord_db[(chr,exon1_stop)]; ji.setLeftExonAnnotations(le); ji.setLeftExonRegionData(exon_region_db[le])
                            re=exon2_coord_db[(chr,exon2_start)]; ji.setRightExonAnnotations(re);  ji.setRightExonRegionData(exon_region_db[re])                                         
                            if le[0] != re[0]: ### Indicates Trans-splicing (e.g., chr7:52,677,568-52,711,750 mouse mm9)
                                ji.setTransSplicing(); #print exon1_stop,le,exon2_start,re,ji.Chr(),ji.Strand()
                            self.both_found+=1; #print 'five',(chr,exon1_stop,exon2_start),exon1_coord_db[(chr,exon1_stop)]
                        else:
                            if (chr,exon1_stop) in exon1_coord_db: ### hence, exon1_stop is known, so report the coordinates of exon2 as novel
                                le=exon1_coord_db[(chr,exon1_stop)]; ji.setLeftExonAnnotations(le)
                                self.one_found+=1; #print 'three',(chr,exon1_stop,exon2_start),exon1_coord_db[(chr,exon1_stop)]
                                novel_exon_coordinates[ji.Chr(),exon2_start] = ji,'left',ji.Exon2Stop() ### Employ this strategy to avoid duplicate exons with differing lengths (mainly an issue if analyzing only exons results)
                                ji.setNovel('side')
                            elif (chr,exon2_start) in exon2_coord_db: ### hence, exon2_start is known, so report the coordinates of exon1 as novel
                                re=exon2_coord_db[(chr,exon2_start)]; ji.setRightExonAnnotations(re) ### In very rare cases, a gene can be assigned here, even though the splice-site is on the opposite strand (not worthwhile filtering out)
                                self.one_found+=1; #print 'three',(chr,exon1_stop,exon2_start),exon1_coord_db[(chr,exon1_stop)]
                                novel_exon_coordinates[ji.Chr(),exon1_stop] = ji,'right',ji.Exon1Start()
                                ji.setNovel('side')
                            else:
                                self.not_found+=1; #if self.not_found < 10: print (chr,exon1_stop,exon2_start)
                                novel_exon_coordinates[ji.Chr(),exon1_stop] = ji,'right',ji.Exon1Start()
                                novel_exon_coordinates[ji.Chr(),exon2_start] = ji,'left',ji.Exon2Stop()
                                ji.setNovel('both')
                            ### We examine reads where one splice-site aligns to a known but the other not, to determine if trans-splicing occurs
                            try: exons_not_identified[chr,ji.Strand()].append((coord,ji))
                            except KeyError: exons_not_identified[chr,ji.Strand()] = [(coord,ji)]
                    """
                    if fl.ExonBedBuildStatus() == 'no':
                        exportNovelJunctions(species,novel_junction_db,condition_count_db,root_dir,dataset_name,'junction') ### Includes known exons
                    """
                    
                    #print self.both_found, ' where both and', self.one_found, 'where one splice-site are known out of',self.both_found+self.one_found+self.not_found
                    #print 'Novel junctions where both splice-sites are known:',self.both_found
                    #print 'Novel junctions where one splice-site is known:',self.one_found
                    #print 'Novel junctions where the splice-sites are not known:',self.not_found
                    clearObjectsFromMemory(exon_region_db); exon_region_db=[] ### Clear memory of this object
        
                    read_aligned_to_gene=0
                    for (chr,strand) in exons_not_identified:
                        if (chr,strand) in self.chr_strand_gene_dbs:
                            chr_gene_locations = self.chr_strand_gene_dbs[chr,strand]
                            chr_reads = exons_not_identified[chr,strand]
                            chr_gene_locations.sort(); chr_reads.sort()
                            ### Set GeneID for each coordinate object (primary and seconardary GeneIDs)
                            read_aligned_to_gene=geneAlign(chr,chr_gene_locations,self.location_gene_db,chr_reads,'no',read_aligned_to_gene)
                    #print read_aligned_to_gene, 'novel junctions aligned to Ensembl genes out of',self.one_found+self.not_found
                    clearObjectsFromMemory(exons_not_identified); exons_not_identified=[] ## Clear memory of this object
                    
                    for key in novel_junction_db:
                        (chr,exon1_stop,exon2_start) = key
                        ji=novel_junction_db[key]
                        if ji.GeneID() == None:
                            try:
                                if ji.SecondaryGeneID() != None:
                                    ### Occurs if mapping is to the 5'UTR of a gene for the left splice-site (novel alternative promoter)
                                    ji.setGeneID(ji.SecondaryGeneID()); ji.setSecondaryGeneID(''); #print key, ji.GeneID(), ji.Strand(), ji.SecondaryGeneID()
                            except Exception: null=[]
                        if ji.GeneID() != None:
                            geneid = ji.GeneID()
                            proceed = 'no'                
                            if ji.SpliceSitesFound() == None: proceed = 'yes'; coordinates = [exon1_stop,exon2_start]
                            elif ji.SpliceSitesFound() == 'left': proceed = 'yes'; coordinates = [exon1_stop,exon2_start]
                            elif ji.SpliceSitesFound() == 'right': proceed = 'yes'; coordinates = [exon1_stop,exon2_start]
                            if proceed == 'yes':
                                for coordinate in coordinates:
                                    if ji.TransSplicing() == 'yes':
                                        #print ji.Chr(),ji.GeneID(), ji.SecondaryGeneID(), ji.Exon1Stop(), ji.Exon2Start()
                                        self.trans_splicing_reads+=1
                                        if ji.checkExonPosition(coordinate) == 'right': geneid = ji.SecondaryGeneID()
                                    exon_data = (coordinate,ji.Chr()+'-'+str(coordinate),'novel')
                                    try: novel_exon_db[geneid].append(exon_data)
                                    except KeyError: novel_exon_db[geneid] = [exon_data]
                                    
                        else: self.junctions_without_exon_gene_alignments+=1
            
                    ### Remove redundant exon entries and store objects    
                    for key in novel_exon_db:
                        exon_data_objects=[]
                        exon_data_list = unique.unique(novel_exon_db[key])
                        exon_data_list.sort()
                        for e in exon_data_list:
                            ed = ExonInfo(e[0],e[1],e[2])
                            exon_data_objects.append(ed)
                        novel_exon_db[key] = exon_data_objects
                        
                    #print self.trans_splicing_reads,'trans-splicing junctions found (two aligning Ensembl genes).'
                    #print self.junctions_without_exon_gene_alignments, 'junctions where neither splice-site aligned to a gene'
                    #if 'X' in search_chr: print len(ens_exon_db),len(ens_exon_db['ENSMUSG00000044424'])
                    alignReadsToExons(novel_exon_db,ens_exon_db,testImport=testImport)
                    
                    ### Link exon annotations up with novel junctions
                    junction_region_db,exons_to_export = annotateNovelJunctions(novel_junction_db,novel_exon_db,exons_to_export)
                    
                    ### Add the exon region data from known Ensembl/UCSC matched junctions to junction_region_db for recipricol junction analysis
                    for key in junction_db:
                        ji=junction_db[key]; jd = ji.ExonAnnotations()
                        try:
                            uid = jd.GeneID()+':'+jd.ExonRegionIDs();  ji.setUniqueID(uid)
                            try: junction_region_db[jd.GeneID()].append((formatID(uid),jd.ExonRegionIDs()))
                            except KeyError: junction_region_db[jd.GeneID()] = [(formatID(uid),jd.ExonRegionIDs())]
                        except AttributeError: null=[] ### Occurs since not all entries in the dictionary are perfect junction matches
                    
                    try: novel_exon_coordinates = collapseNoveExonBoundaries(novel_exon_coordinates,root_dir+dataset_name) ### Joins inferred novel exon-IDs (5' and 3' splice sites) from adjacent and close junction predictions
                    except Exception: pass ### No errors encountered before
                    if self.exonBedBuildStatus == 'yes':
                        ### Append to the exported BED format exon coordinate file
                        bedfile = exportNovelExonToBedCoordinates(species,novel_exon_coordinates,chr_status,searchChr=searchchr)
                    
                    ### Identify reciprocol junctions and retrieve splice-event annotations for exons and inclusion junctions
                    
                    junction_annotations,critical_exon_annotations = JunctionArray.inferJunctionComps(species,('RNASeq',junction_region_db,root_dir),searchChr=searchchr)
                    clearObjectsFromMemory(junction_region_db); junction_region_db=[]   
                    
                    ### Reformat these dictionaries to combine annotations from multiple reciprocol junctions
                    junction_annotations = combineExonAnnotations(junction_annotations)
                    critical_exon_annotations = combineExonAnnotations(critical_exon_annotations)
                    
                if 'exon' in biotypes:
                    if testImport == 'yes':
                        print len(unmapped_exon_db),'exon genomic locations imported.'
                    ### Create a new dictionary keyed by chromosome and strand
                    exons_not_aligned={}
                    for (chr,exon1_stop,exon2_start) in unmapped_exon_db:
                        ji = unmapped_exon_db[(chr,exon1_stop,exon2_start)]
                        coord = [exon1_stop,exon2_start]; coord.sort()
                        try: exons_not_aligned[chr,ji.Strand()].append((coord,ji))
                        except KeyError: exons_not_aligned[chr,ji.Strand()] = [(coord,ji)]
                                                          
                    read_aligned_to_gene=0
                    for (chr,strand) in exons_not_aligned:
                        if (chr,strand) in self.chr_strand_gene_dbs:
                            chr_gene_locations = self.chr_strand_gene_dbs[chr,strand]
                            chr_reads = exons_not_aligned[chr,strand]
                            chr_gene_locations.sort(); chr_reads.sort()
                            read_aligned_to_gene=geneAlign(chr,chr_gene_locations,self.location_gene_db,chr_reads,'no',read_aligned_to_gene)
                            
                    #print read_aligned_to_gene, 'exons aligned to Ensembl genes out of',self.one_found+self.not_found
                    
                    align_exon_db={}; exons_without_gene_alignments={}; multigene_exon=0
                    for key in unmapped_exon_db:
                        (chr,exon1_stop,exon2_start) = key
                        ji=unmapped_exon_db[key]
                        if ji.GeneID() == None:
                            try:
                                if ji.SecondaryGeneID() != None:
                                    ### Occurs if mapping outside known exon boundaries for one side of the exon
                                    ji.setGeneID(ji.SecondaryGeneID()); ji.setSecondaryGeneID(''); #print key, ji.GeneID(), ji.Strand(), ji.SecondaryGeneID()
                            except Exception: null=[]
                        else:
                            if 'ENS' in ji.JunctionID():
                                if ji.GeneID() not in ji.JunctionID(): ### Hence, there were probably two overlapping Ensembl genes and the wrong was assigned based on the initial annotations
                                    original_geneid = string.split(ji.JunctionID(),':')[0]
                                    if original_geneid in ens_exon_db: ji.setGeneID(original_geneid) #check if in ens_exon_db (since chromosome specific)
                                    
                        if ji.GeneID() != None:
                            geneid = ji.GeneID()
                            coordinates = [exon1_stop,exon2_start]
                            for coordinate in coordinates:
                                if ji.TransSplicing() != 'yes': ### This shouldn't occur for exons
                                    exon_data = (coordinate,ji.Chr()+'-'+str(coordinate),'novel')
                                    try: align_exon_db[geneid].append(exon_data)
                                    except KeyError: align_exon_db[geneid] = [exon_data]
                                else:
                                    multigene_exon+=1 ### Shouldn't occur due to a fix in the gene-alignment method which will find the correct gene on the 2nd interation
                        else: exons_without_gene_alignments[key]=ji; self.exons_without_gene_alignment_count+=1
        
                    ### Remove redundant exon entries and store objects (this step may be unnecessary)
                    for key in align_exon_db:
                        exon_data_objects=[]
                        exon_data_list = unique.unique(align_exon_db[key])
                        exon_data_list.sort()
                        for e in exon_data_list:
                            ed = ExonInfo(e[0],e[1],e[2])
                            exon_data_objects.append(ed)
                        align_exon_db[key] = exon_data_objects
                        
                    #print self.exons_without_gene_alignment_count, 'exons where neither aligned to a gene'
                    #if self.exons_without_gene_alignment_count>3000: print 'NOTE: Poor mapping of these exons may be due to an older build of\nEnsembl than the current version. Update BAMtoBED mappings to correct.'
            
                    begin_time = time.time()
                    alignReadsToExons(align_exon_db,ens_exon_db)
                    
                    end_time = time.time()
                    if testImport == 'yes':
                        print 'Exon sequences aligned to exon regions in',int(end_time-begin_time),'seconds'
        
                    ### Combine the start and end region alignments into a single exon annotation entry
                    combineDetectedExons(unmapped_exon_db,align_exon_db,novel_exon_db)
                    clearObjectsFromMemory(unmapped_exon_db); clearObjectsFromMemory(align_exon_db); clearObjectsFromMemory(novel_exon_db)
                    unmapped_exon_db=[]; align_exon_db=[]; novel_exon_db=[]
                    """            
                    if fl.ExonBedBuildStatus() == 'no':
                        exportNovelJunctions(species,exons_without_gene_alignments,condition_count_db,root_dir,dataset_name,'exon') ### Includes known exons
                    """
                    clearObjectsFromMemory(exons_without_gene_alignments); exons_without_gene_alignments=[]
                    
                ### Export both exon and junction annotations
                
                if 'junction' in biotypes:
                    ### Export the novel user exon annotations    
                    exportDatasetLinkedExons(species,exons_to_export,critical_exon_annotations,root_dir,testImport=testImport,searchChr=searchchr)            
                            
                ### Export the novel user exon-junction annotations (original junction_db objects updated by above processing)
                exportDatasetLinkedJunctions(species,junction_db,junction_annotations,root_dir,testImport=testImport,searchChr=searchchr)
        
                ### Clear memory once results are exported (don't want to delete actively used objects)
                if 'junction' in biotypes:
                    clearObjectsFromMemory(exons_to_export); clearObjectsFromMemory(critical_exon_annotations)
                    clearObjectsFromMemory(novel_junction_db); novel_junction_db=[]
                    clearObjectsFromMemory(novel_exon_coordinates); novel_exon_coordinates=[]    
                    exons_to_export=[]; critical_exon_annotations=[]
                    clearObjectsFromMemory(exon1_coord_db); clearObjectsFromMemory(exon2_coord_db)
                    exon1_coord_db=[]; exon2_coord_db=[]
                if 'exon' in biotypes:            
                    clearObjectsFromMemory(exons_not_aligned); exons_not_aligned=[]
                    clearObjectsFromMemory(ens_exon_db); ens_exon_db=[]
                                
                ### Add chromsome specific junction_db data to a simple whole genome dictionary
                for key in junction_db:
                    ji = junction_db[key]
                    if ji.GeneID()!=None and ji.UniqueID()!=None: self.junction_simple_db[key]=ji.UniqueID()
        
                #returnLargeGlobalVars()
                clearObjectsFromMemory(junction_db); clearObjectsFromMemory(junction_annotations)
                junction_db=[]; junction_annotations=[]; chr_reads=[]
    
                for biotype in biotypes:
                    ### Import Read Counts (do this last to conserve memory)
            
                    if platformType == 'RNASeq':
    
                         condition_count_db,exon_len_db,biotypes2,algorithms = importBEDFile(bed_dir,root_dir,species,normalize_feature_exp,getReads=True,searchChr=searchchr,getBiotype=biotype,testImport=testImport)
                    else:
                         condition_count_db,exon_len_db,biotypes2,algorithms = importExpressionMatrix(bed_dir,root_dir,species,fl,'yes',getBiotype=biotype)
                    ###First export original counts, rather than quantile normalized or RPKM
                    self.exportJunctionCounts(species,self.junction_simple_db,exon_len_db,condition_count_db,root_dir,dataset_name,biotype,'counts',searchChr=searchchr)
                    clearObjectsFromMemory(condition_count_db); clearObjectsFromMemory(exon_len_db); condition_count_db=[]; exon_len_db=[]
            if analysisMode == 'commandline':
                    print 'finished parsing data for chromosome:',search_chr ### Unix platforms are not displaying the progress in real-time
            else:
                pass #print "*",
            try: queue.put([self.biotypes_store, self.known_count, self.novel_junction_count, self.trans_splicing_reads, self.junctions_without_exon_gene_alignments, self.exons_without_gene_alignment_count])
            except Exception:
                ### If queue is not a multiprocessing object
                queue = [self.biotypes_store, self.known_count, self.novel_junction_count, self.trans_splicing_reads, self.junctions_without_exon_gene_alignments, self.exons_without_gene_alignment_count]
            return queue
        
        except Exception:
            print traceback.format_exc()
            error(traceback.format_exc())
            multiprocessing.log_to_stderr().setLevel(logging.DEBUG)
            raise
    
    def exportJunctionCounts(self,species,junction_simple_db,exon_len_db,condition_count_db,root_dir,dataset_name,biotype,count_type,searchChr=None):
        if 'exp.' not in dataset_name: dataset_name = 'exp.'+dataset_name
        if '.txt' not in dataset_name: dataset_name+='.txt'
        export_path = root_dir+'ExpressionInput/'+dataset_name
        if count_type == 'counts':
            export_path = string.replace(export_path,'exp.','counts.') ### separately export counts
            if searchChr !=None:
                export_path = string.replace(export_path,'ExpressionInput','ExpressionInput/Counts')
                export_path = string.replace(export_path,'.txt','.'+searchChr+'.txt')

        self.countsFile = export_path
        
        if self.testImport == 'yes':
            print 'Writing',export_path
        export_data,status = AppendOrWrite(export_path)
        
        if status == 'not found':
            title = ['AltAnalyze_ID']
            for condition in condition_count_db: title.append(condition)
            export_data.write(string.join(title,'\t')+'\n')
        
        for key in self.junction_simple_db:
            chr,exon1_stop,exon2_start = key
            if biotype == 'junction':
                coordinates = chr+':'+str(exon1_stop)+'-'+str(exon2_start)
            elif biotype == 'exon':
                coordinates = chr+':'+str(exon1_stop-1)+'-'+str(exon2_start+1)
            try:
                null=exon_len_db[key]
                if count_type == 'counts': values = [self.junction_simple_db[key]+'='+coordinates]
                else: values = [self.junction_simple_db[key]]
                for condition in condition_count_db: ###Memory crash here
                    count_db = condition_count_db[condition]
                    try: read_count = count_db[key]
                    except KeyError: read_count = '0'
                    values.append(read_count)
                export_data.write(string.join(values,'\t')+'\n')
            except Exception: null=[]
        export_data.close()
    
    def countsDir(self):
        return self.countsFile
    
def fastRPKMCalculate(counts_file):
    export_path = string.replace(counts_file,'counts.','exp.')
    export_data = export.ExportFile(export_path) ### Write this new file
    fn=filepath(counts_file); header=True
    exon_sum_array=[]; junction_sum_array=[]
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if header:
            samples = t[1:]
            header=False
            exon_sum_array=[0]*len(samples)
            junction_sum_array=[0]*len(samples)
        else:
            try: values = map(float,t[1:])
            except Exception: print t; sys.exit()
            ### get the total reads/sample
            if '-' in string.split(t[0],'=')[0]:
                junction_sum_array = [sum(value) for value in zip(*[junction_sum_array,values])]
            else:
                exon_sum_array = [sum(value) for value in zip(*[exon_sum_array,values])]
                
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",category=RuntimeWarning) ### hides warnings associated with Scipy for n=1 sample comparisons
        jatr=Average(junction_sum_array) # Average of the total maped reads
        eatr=Average(exon_sum_array) # Average of the total maped reads
    
    header=True
    c=math.pow(10.0,9.0)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if header:
            export_data.write(line) ### Write header
            header=False
        else:
            exon_id,coordinates = string.split(t[0],'=')
            coordinates = string.split(coordinates,':')[1]
            coordinates = string.split(coordinates,'-')
            l=abs(int(coordinates[1])-int(coordinates[0])) ### read-length
                        
            read_counts = map(lambda x: int(x)+1, t[1:])
            if '-' in exon_id:
                count_stats = zip(read_counts,junction_sum_array)
                atr = jatr
                l=60
            else:
                count_stats = zip(read_counts,exon_sum_array)
                atr = eatr
            values=[]
            
            #rpkm = map(lambda (r,t): c*(r/(t*l)), count_stats) ### Efficent way to convert to rpkm, but doesn't work for 0 counts
            for (r,t) in count_stats:
                if r == 1: ###This adjustment allows us to obtain more realist folds where 0 is compared and use log2
                    t = atr
                try:
                    rpkm = str(c*(r/(t*l)))
                    #print c,r,t,l,exon_id,rpkm;sys.exit()
                    values.append(rpkm)
                except Exception,e:
                    print e
                    print t[0]
                    print 'Error Encountered... Exon or Junction of zero length encoutered... RPKM failed... Exiting AltAnalyze.'
                    print 'This error may be due to inconsistent file naming. If both exon and junction sample data is present, make sure they are named propperly.'
                    print 'For example: cancer1__exon.bed, cancer1__junction.bed (double underscore required to match these samples up)!'
                    print [r,t,l];k=1; forceError
            values = string.join([exon_id]+values,'\t')+'\n'
            export_data.write(values)
    export_data.close()
            
def mergeCountFiles(counts_file1,counts_file2):
    ### Used internally to merge count files that are very large and too time-consuming to recreate (regenerate them)
    export_path = string.replace(counts_file2,'counts.','temp-counts.')
    export_data = export.ExportFile(export_path) ### Write this new file
    fn=filepath(counts_file1); header=True
    count_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if header:
            samples = t[1:]
            header=False
            si = samples.index('H9.102.2.5.bed')+1
            
        else:
            try: value = t[si]
            except Exception: print t; sys.exit()
            ### get the total reads/sample
            count_db[t[0]] = value

    fn=filepath(counts_file2); header=True
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if header:
            samples = t[1:]
            header=False
            si = samples.index('H9.102.2.5.bed')+1
            export_data.write(line)
        else:
            try: t[si] = count_db[t[0]]
            except Exception: pass ### keep the current value
            export_data.write(string.join(t,'\t')+'\n')
    export_data.close()
    
def importRawCountData(filename,expressed_gene_exon_db):
    """ Identifies exons or junctions to evaluate gene-level expression. This function, as it is currently written:
    1) examines the RPKM and original read counts associated with all exons
    2) removes exons/junctions that do not meet their respective RPKM AND read count cutoffs
    3) returns ONLY those exons and genes deemed expressed, whether constitutive selected or all exons
    """
    
    ### Get expression values for exon/junctions to analyze
    seq_ids_to_import={}
    for gene in expressed_gene_exon_db:
        for exonid in expressed_gene_exon_db[gene]: seq_ids_to_import[exonid]=[]
            
    ### Define thresholds
    exon_exp_threshold = UserOptions.ExonExpThreshold()
    junction_exp_threshold = UserOptions.JunctionExpThreshold()
    exon_rpkm_threshold = UserOptions.ExonRPKMThreshold()

    gene_rpkm_threshold = UserOptions.RPKMThreshold()
    gene_exp_threshold = UserOptions.GeneExpThreshold()
    
    ### Import RPKM normalized expression values               
    fn=filepath(filename); x=0; rpkm_dbase={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: array_names = t[1:]; x=1
        else:
            exon_id=t[0]
            max_count=max(map(float,t[1:]))
            if max_count>=exon_rpkm_threshold: rpkm_dbase[exon_id]=[] ### Only retain exons/junctions meeting the RPKM threshold
            
    ### Import non-normalized original counts                
    counts_filename = string.replace(filename,'exp.','counts.')
    fn=filepath(counts_filename); x=0; exp_dbase={}
    all_exp_features={} ### Don't filter for only gene-expression reporting
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: array_names = t[1:]; x=1
        else:
            exon_id,coordinates = string.split(t[0],'=')
            coordinates = string.split(coordinates,':')[1]
            coordinates = string.split(coordinates,'-')
            length=abs(int(coordinates[1])-int(coordinates[0]))
            max_count=max(map(float,t[1:])); proceed = 'no'
            if '-' in exon_id:
                length = 60.0
                if max_count>=junction_exp_threshold:
                    ### Only considered when exon data is not present in the analysis
                    proceed = 'yes'
            elif max_count>=exon_exp_threshold: proceed = 'yes'
            if proceed == 'yes' and exon_id in rpkm_dbase: ### Ensures that the maximum sample (not group) user defined count threshold is achieved at the exon or junction-level
                all_exp_features[exon_id]=None
                if exon_id in seq_ids_to_import:### Forces an error if not in the steady-state pre-determined set (CS or all-exons) - INCLUDE HERE TO FILTER ALL FEATURES
                    exp_dbase[exon_id] = t[1:],length ### Include sequence length for normalization

    for exon in exp_dbase: array_count = len(exp_dbase[exon][0]); break
    try:null=array_count
    except Exception:
        print 'No exons or junctions considered expressed (based user thresholds). Exiting analysis.'; force_exit
    return exp_dbase, all_exp_features, array_count

def importNormalizedCountData(filename,expressed_gene_exon_db):
    ### Get expression values for exon/junctions to analyze
    seq_ids_to_import={}
    for gene in expressed_gene_exon_db:
        for exonid in expressed_gene_exon_db[gene]: seq_ids_to_import[exonid]=[]
            
    ### Define thresholds
    exon_exp_threshold = UserOptions.ExonExpThreshold()
    junction_exp_threshold = UserOptions.JunctionExpThreshold()
    exon_rpkm_threshold = UserOptions.ExonRPKMThreshold()
    
    gene_rpkm_threshold = UserOptions.RPKMThreshold()
    gene_exp_threshold = UserOptions.GeneExpThreshold()
    
    ### Import non-normalized original counts                
    fn=filepath(filename); x=0; exp_dbase={}
    all_exp_features={} ### Don't filter for only gene-expression reporting
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0: array_names = t[1:]; x=1
        else:
            exon_id=t[0]; proceed = 'no'
            max_count=max(map(float,t[1:]))
            if '-' in exon_id:
                if max_count>=junction_exp_threshold: proceed = 'yes'
            elif max_count>=exon_exp_threshold: proceed = 'yes'
            if proceed == 'yes': ### Ensures that the maximum sample (not group) user defined count threshold is achieved at the exon or junction-level
                all_exp_features[exon_id]=None
                if exon_id in seq_ids_to_import: ### If a "constitutive" or exon-level feature (filter missing prior to 2.0.8 - bug)
                    exp_dbase[exon_id] = t[1:],0  ### Add the zero just to comply with the raw count input format (indicates exon length)

    for exon in exp_dbase: array_count = len(exp_dbase[exon][0]); break        
    return exp_dbase, all_exp_features, array_count

def obtainGeneCounts(expressed_gene_exon_db,exp_dbase,array_count,normalize_feature_exp):
    ###Calculate avg expression for each sample for each exon (using constitutive or all exon values)

    steady_state_db={}
    for gene in expressed_gene_exon_db:
        x = 0; gene_sum=0
        exon_list = expressed_gene_exon_db[gene]
        while x < array_count:
            exp_list=[]; len_list=[]
            for exon in exon_list:
                try:
                    exp_val = exp_dbase[exon][0][x]
                    if normalize_feature_exp == 'RPKM':
                        ### Decided to include all exons, expressed or not to prevent including lowly expressed exons that are long, that can bias the expression call
                        #if float(exp_val) != 0: ### Here, we use the original raw count data, whereas above is the adjusted quantile or raw count data
                        exp_list.append(exp_val); len_list.append(exp_dbase[exon][1]) ### This is for RNASeq -> don't include undetected exons - made in v.204
                    else: exp_list.append(exp_val) #elif float(exp_val) != 1:
                except KeyError: null =[] ###occurs if the expression exon list is missing some of these exons
            try:
                if len(exp_list)==0:
                    for exon in exon_list:
                        try:
                            exp_list.append(exp_dbase[exon][0][x]); len_list.append(exp_dbase[exon][1])
                            #kill
                        except KeyError: null=[] ### Gene entries will cause this error, since they are in the database but not in the count file
                if normalize_feature_exp == 'RPKM':
                    sum_const_exp=sum(map(float,exp_list)); gene_sum+=sum_const_exp
                    sum_length=sum(len_list) ### can have different lengths for each sample, since only expressed exons are considered
                    ### Add only one avg-expression value for each array, this loop
                    try: steady_state_db[gene].append((sum_const_exp,sum_length))
                    except KeyError: steady_state_db[gene] = [(sum_const_exp,sum_length)]
                else:
                    avg_const_exp=Average(exp_list)
                    if avg_const_exp != 1: gene_sum+=avg_const_exp
                    ### Add only one avg-expression value for each array, this loop
                    try: steady_state_db[gene].append(avg_const_exp)
                    except KeyError: steady_state_db[gene] = [avg_const_exp]
            except Exception: null=[] ### Occurs when processing a truncated dataset (for testing usually) - no values for the gene should be included
            x += 1
        if gene_sum==0:
            try:
                del steady_state_db[gene] ### Hence, no genes showed evidence of expression (most critical for RNA-Seq)
            except Exception: null=[] ### Error occurs when a gene is added to the database from self.location_gene_db, but is not expressed
            
    return steady_state_db

def returnLargeGlobalVars():
    ### Prints all large global variables retained in memory (taking up space)
    all = [var for var in globals() if (var[:2], var[-2:]) != ("__", "__")]
    for var in all:
        try:
            if len(globals()[var])>1:
                print var, len(globals()[var])
        except Exception: null=[]
        
def clearObjectsFromMemory(db_to_clear):
    db_keys={}
    for key in db_to_clear: db_keys[key]=[]
    for key in db_keys:
        try: del db_to_clear[key]
        except Exception: 
            try:
                for i in key: del i ### For lists of tuples
            except Exception: del key ### For plain lists

def verifyFile(filename):
    status = 'not found'
    try:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines(): status = 'found';break
    except Exception: status = 'not found'
    return status

def AppendOrWrite(export_path):
    export_path = filepath(export_path)
    status = verifyFile(export_path)
    if status == 'not found':
        export_data = export.ExportFile(export_path) ### Write this new file
    else:
        export_data = open(export_path,'a') ### Appends to existing file
        
    return export_data, status

def quantileNormalizationSimple(condition_count_db):
    ### Basic quantile normalization method (average ranked expression values)

    ### Get all junction or exon entries
    key_db={}
    for condition in condition_count_db:
        count_db = condition_count_db[condition]
        for key in count_db: key_db[key]=[]
            
    condition_unnormalized_db={}
    for key in key_db:
        ### Only look at the specific biotype of interest for each normalization
        for condition in condition_count_db:
            count_db = condition_count_db[condition]
            try:
                count = float(count_db[key])+1 ###This adjustment allows us to obtain more realist folds where 0 is compared and use log2
                count_db[key] = [] ### Set equal to null as a temporary measure to save memory
            except KeyError: count = 1.00 ###Was zero, but needs to be one for more realistic log2 fold calculations                
            ### store the minimal information to recover the original count and ID data prior to quantile normalization
            try: condition_unnormalized_db[condition].append([count,key])
            except Exception: condition_unnormalized_db[condition]=[[count,key]]
                    
    quantile_normalize_db={}; key_db={}
    for condition in condition_unnormalized_db:
        condition_unnormalized_db[condition].sort() ### Sort lists by count number
        rank=0 ### thus, the ID is the rank order of counts
        for (count,key) in condition_unnormalized_db[condition]:
            try: quantile_normalize_db[rank].append(count)
            except KeyError: quantile_normalize_db[rank] = [count]
            rank+=1
            
    ### Get the average value for each index
    for rank in quantile_normalize_db:
        quantile_normalize_db[rank] = Average(quantile_normalize_db[rank])

    for condition in condition_unnormalized_db:
        rank=0
        count_db = condition_count_db[condition]
        for (count,key) in condition_unnormalized_db[condition]:
            avg_count = quantile_normalize_db[rank]
            rank+=1
            count_db[key] = str(avg_count) ### re-set this value to the normalized value
    
    try:
        clearObjectsFromMemory(condition_unnormalized_db); condition_unnormalized_db =  []
        clearObjectsFromMemory(quantile_normalize_db); quantile_normalize_db = []
    except Exception: None
    
    return condition_count_db

def combineExonAnnotations(db):
    for i in db:
        list1=[]; list2=[]
        for (junctions,splice_event) in db[i]:
            list1.append(junctions); list2.append(splice_event)
        junctions = EnsemblImport.combineAnnotations(list1)
        splice_event = EnsemblImport.combineAnnotations(list2)
        db[i] = junctions,splice_event
    return db

def formatID(id):
    ### JunctionArray methods handle IDs with ":" different than those that lack this
    return string.replace(id,':','@')

            
def getChromosomeStrandCoordinates(species,testImport):
    ### For novel junctions with no known-splice site, map to genes
    gene_location_db = EnsemblImport.getEnsemblGeneLocations(species,'RNASeq','key_by_array')
    
    chr_strand_gene_db = {}; location_gene_db = {}; chromosome_names={}; all_chromosomes={}
    for gene in gene_location_db:
        chr,strand,start,end = gene_location_db[gene]
        location_gene_db[chr,int(start),int(end)] = gene,strand
        try: chr_strand_gene_db[chr,strand].append((int(start),int(end)))
        except KeyError: chr_strand_gene_db[chr,strand] = [(int(start),int(end))]
        if testImport == 'yes':
            if chr=='chr1': chromosome_names[chr]=[]
            if chr=='chr19': chromosome_names[chr]=[] ### Gene rich chromosome
            if chr=='chrMT': chromosome_names[chr]=[] ### Gene rich chromosome
        elif len(chr)<7: chromosome_names[chr]=[]
        all_chromosomes[chr]=[]
    ### Some organisms aren't organized into classical chromosomes (why I don't know)
    if len(chromosome_names)<10 and len(all_chromosomes)>9 and testImport=='no': chromosome_names = all_chromosomes
    return chr_strand_gene_db,location_gene_db,chromosome_names,gene_location_db

def exportDatasetLinkedExons(species,exons_to_export,critical_exon_annotations,root_dir,testImport=None,searchChr=None):
    export_path = root_dir+'AltDatabase/'+species+'/RNASeq/'+species + '_Ensembl_exons.txt'
    if searchChr != None:
        export_path = string.replace(export_path,'RNASeq/'+species,'RNASeq/exons/'+species)
        export_path = string.replace(export_path,'.txt','.'+searchChr+'.txt')
    
    if testImport == 'yes': print 'Writing',export_path
    export_data,status = AppendOrWrite(export_path)
    
    if status == 'not found': 
        export_title = ['AltAnalyzeID','exon_id','ensembl_gene_id','transcript_cluster_id','chromosome','strand','probeset_start','probeset_stop']
        export_title +=['class','constitutive_probeset','ens_exon_ids','ens_constitutive_status','exon_region','exon-region-start(s)','exon-region-stop(s)','splice_events','splice_junctions']
        export_title = string.join(export_title,'\t')+'\n'; export_data.write(export_title)

    ### We stored these in a dictionary to make sure each exon is written only once and so we can organize by gene
    exons_to_export_list=[]
    for key in exons_to_export:
        ed = exons_to_export[key]
        exons_to_export_list.append((key,ed))
    exons_to_export_list.sort()
    
    for (key,ed) in exons_to_export_list:
        constitutive_call = 'no'; ens_constitutive_status = '0'
        try:
            red = ed.ExonRegionData()
            exon_region = ed.ExonRegionID()
            start = str(ed.ReadStart()); stop = start
            if '-' not in exon_region and '_' not in exon_region: annotation = 'known'
            else: annotation = 'novel'      
        except Exception:
            red = ed ### For annotated exons, no difference in the annotations
            exon_region = ed.ExonRegionIDs()
            start = str(red.ExonStart()); stop = str(red.ExonStop())
            constitutive_call = red.Constitutive()
            if constitutive_call == 'yes': ens_constitutive_status = '1'
            annotation = 'known'            
        uid = red.GeneID()+':'+exon_region
        splice_events = red.AssociatedSplicingEvent(); splice_junctions = red.AssociatedSplicingJunctions()
        if uid in critical_exon_annotations:
            splice_junctions,splice_events = critical_exon_annotations[uid]
        export_values = [uid, exon_region, red.GeneID(), '', red.Chr(), red.Strand(), start, stop, annotation, constitutive_call, red.ExonID(), ens_constitutive_status]
        export_values+= [exon_region, str(red.ExonStart()), str(red.ExonStop()), splice_events, splice_junctions]
        export_values = string.join(export_values,'\t')+'\n'; export_data.write(export_values)
    export_data.close()
    
def exportNovelJunctions(species,novel_junction_db,condition_count_db,root_dir,dataset_name,biotype):
    
    if 'exp.' not in dataset_name: dataset_name = 'exp.'+dataset_name
    if '.txt' not in dataset_name: dataset_name+='.txt'
    dataset_name = string.replace(dataset_name,'exp','novel')
    dataset_name = string.replace(dataset_name,'.txt','.'+biotype+'.txt')
    
    export_path = root_dir+'ExpressionInput/'+dataset_name
    export_data,status = AppendOrWrite(export_path)
    
    if status == 'not found': 
        title = ['chr','strand','start','stop','start Ensembl','end Ensembl','known start', 'known end']
        for condition in condition_count_db: title.append(condition)
        export_data.write(string.join(title,'\t')+'\n')
    
    for key in novel_junction_db:
        ji = novel_junction_db[key]
        try: gene1 = str(ji.GeneID())
        except Exception: gene1=''
        try: gene2 = str(ji.SecondaryGeneID())
        except Exception: gene2 = 'None'
        try: le = str(ji.LeftExonAnnotations())
        except Exception: le = ''
        try: re = str(ji.RightExonAnnotations())
        except Exception: re = ''
        if biotype == 'junction':
            values = [ji.Chr(), ji.Strand(), str(ji.Exon1Stop()), str(ji.Exon2Start())]
        elif biotype == 'exon':
            values = [ji.Chr(), ji.Strand(), str(ji.Exon1Stop()-1), str(ji.Exon2Start()+1)] ### correct for initial adjustment
        values += [gene1,gene2,le,re]
        for condition in condition_count_db:
            count_db = condition_count_db[condition]
            try: read_count = count_db[key]
            except KeyError: read_count = '0'
            values.append(read_count)
        export_data.write(string.join(values,'\t')+'\n')
    export_data.close()

def exportDatasetLinkedGenes(species,gene_location_db,root_dir):
    """Include an entry for gene IDs to include constitutive expression for RPKM normalized data"""

    export_path = root_dir+'AltDatabase/'+species+'/RNASeq/'+species + '_Ensembl_junctions.txt'
    export_data,status = AppendOrWrite(export_path)
    for gene in gene_location_db:
        chr,strand,start,end = gene_location_db[gene]
        export_values = [gene, 'E0.1',gene, '', chr, strand, str(start), str(end), 'known', 'yes', gene, '1']
        export_values+= ['E0.1', str(start), str(end), '', '']
        export_values = string.join(export_values,'\t')+'\n'; export_data.write(export_values)
    export_data.close()
               
def exportDatasetLinkedJunctions(species,junction_db,junction_annotations,root_dir,testImport=False,searchChr=None):
    export_path = root_dir+'AltDatabase/'+species+'/RNASeq/'+species + '_Ensembl_junctions.txt'
    if searchChr != None:
        
        export_path = string.replace(export_path,'RNASeq/'+species,'RNASeq/junctions/'+species)
        export_path = string.replace(export_path,'.txt','.'+searchChr+'.txt')
        
    if testImport == 'yes': print 'Writing',export_path
    export_data,status = AppendOrWrite(export_path)

    if status == 'not found': 
        export_title = ['AltAnalyzeID','exon_id','ensembl_gene_id','transcript_cluster_id','chromosome','strand','probeset_start','probeset_stop']
        export_title +=['class','constitutive_probeset','ens_exon_ids','ens_constitutive_status','exon_region','exon-region-start(s)','exon-region-stop(s)','splice_events','splice_junctions']
        export_title = string.join(export_title,'\t')+'\n'; export_data.write(export_title)

    for key in junction_db:
        (chr,exon1_stop,exon2_start) = key
        ji=junction_db[key]
        #print key, ji.UniqueID(), ji.GeneID()
        if ji.GeneID()!=None and ji.UniqueID()!=None:
            if ji.UniqueID() in junction_annotations: ### Obtained from JunctionArray.inferJunctionComps()
                junctions,splice_events = junction_annotations[ji.UniqueID()]
                if ji.TransSplicing() == 'yes':
                    if len(splice_events)>0: splice_events+= '|trans-splicing'
                    else: splice_events = 'trans-splicing'
                ji.setAssociatedSplicingEvent(splice_events); ji.setAssociatedSplicingJunctions(junctions)
            elif ji.TransSplicing() == 'yes':
                ji.setAssociatedSplicingEvent('trans-splicing')
            try:
                try: constitutive_call = ji.Constitutive()
                except Exception:
                    jd = ji.ExonAnnotations()
                    constitutive_call = jd.Constitutive()
                if constitutive_call == 'yes': ens_constitutive_status = '1'
                else: ens_constitutive_status = '0'
                annotation = 'known'
            except Exception:
                constitutive_call = 'no'; ens_constitutive_status = '0'; annotation = 'novel'
            if 'I' in ji.ExonRegionID() or 'U' in ji.ExonRegionID() or '_' in ji.ExonRegionID():
                annotation = 'novel' ### Not previously indicated well (as I remember) for exon-level reads - so do this
            export_values = [ji.UniqueID(), ji.ExonRegionID(), ji.GeneID(), '', ji.Chr(), ji.Strand(), str(ji.Exon1Stop()), str(ji.Exon2Start()), annotation, constitutive_call, ji.ExonID(), ens_constitutive_status]
            export_values+= [ji.ExonRegionID(), str(ji.Exon1Stop()), str(ji.Exon2Start()), ji.AssociatedSplicingEvent(), ji.AssociatedSplicingJunctions()]
            export_values = string.join(export_values,'\t')+'\n'; export_data.write(export_values)
    export_data.close()
       
def combineDetectedExons(unmapped_exon_db,align_exon_db,novel_exon_db):
    ### Used for exon alignments (both start position and end position aligned to exon/intron/UTR regions)
    ### Reformat align_exon_db to easily lookup exon data
    aligned_exon_lookup_db={}
    for gene in align_exon_db:
        for ed in align_exon_db[gene]:
            aligned_exon_lookup_db[gene,ed.ReadStart()]=ed
            #if gene == 'ENSMUSG00000064181': print ed.ReadStart(),ed.ExonRegionID()
            
    ### Reformat novel_exon_db to easily lookup exon data - created from junction analysis (rename above exons to match novel junctions)
    novel_exon_lookup_db={}
    for gene in novel_exon_db:
        for ed in novel_exon_db[gene]:
            try:
                ### Only store exons that are found in the novel exon file
                null = aligned_exon_lookup_db[gene,ed.ReadStart()+1] ### offset introduced on import
                novel_exon_lookup_db[gene,ed.ReadStart()+1]=ed
            except Exception: null=[]
            try:
                ### Only store exons that are found in the novel exon file
                null = aligned_exon_lookup_db[gene,ed.ReadStart()-1] ### offset introduced on import
                novel_exon_lookup_db[gene,ed.ReadStart()-1]=ed
            except Exception: null=[]
                    
    ### Lookup the propper exon region ID and gene ID to format the unique ID and export coordinates
    x = 0
    for key in unmapped_exon_db:
        (chr,exon1_stop,exon2_start) = key
        ji=unmapped_exon_db[key]
        proceed = 'no'
        if ji.GeneID() != None:
            e1 = (ji.GeneID(),exon1_stop)
            e2 = (ji.GeneID(),exon2_start)
            exon_info=[]; override_annotation = None; found=[]
            try: null = aligned_exon_lookup_db[e1]; found.append(1)
            except Exception: null=[]
            try: null = aligned_exon_lookup_db[e2]; found.append(2)
            except Exception: null=[]
            try: null = novel_exon_lookup_db[e1]; override_annotation = 1
            except Exception:
                try: null = novel_exon_lookup_db[e2]; override_annotation = 2
                except Exception: null=[]          
            if len(found)>0:
                ### Below is not the simplist way to do this, but should be the fastest
                if 1 in found: exon_info.append(aligned_exon_lookup_db[e1])
                if 2 in found: exon_info.append(aligned_exon_lookup_db[e2])
                if len(exon_info) == 2: ed1,ed2 = exon_info
                else:
                    ed1 = exon_info[0]; ed2 = ed1; x+=1 ### if only one splice site aligned to a gene region (shouldn't occur)
                    if x == 2: null=[]; #print 'SOME EXONS FOUND WITH ONLY ONE ALIGNING POSITION...',key,ji.GeneID(),ed1.ExonRegionID(),e1,e2
                try: red1 = ed1.ExonRegionData(); red2 = ed2.ExonRegionData()
                except Exception:
                    print [ji.GeneID(), ji.Chr(), key]
                    print e1, e2
                    try: print ed1.ExonRegionData()
                    except Exception: 'ed1 failed'
                    try: print ed2.ExonRegionData()
                    except Exception: 'ed2 failed'
                    kill
                region1 = ed1.ExonRegionID(); region2 = ed2.ExonRegionID()
                #print region1,region2,ji.GeneID(),ji.Chr(),ji.Strand()
             
                try: splice_junctions = EnsemblImport.combineAnnotations([red1.AssociatedSplicingJunctions(),red2.AssociatedSplicingJunctions()])
                except Exception: print red1, red2;sys.exit()
                splice_events = EnsemblImport.combineAnnotations([red1.AssociatedSplicingEvent(),red2.AssociatedSplicingEvent()])
                ji.setAssociatedSplicingJunctions(splice_junctions)
                ji.setAssociatedSplicingEvent(splice_events)
                ens_exon_ids = EnsemblImport.combineAnnotations([red1.ExonID(),red2.ExonID()])
                ji.setExonID(ens_exon_ids)
                if red1.Constitutive() == 'yes' or red2.Constitutive() == 'yes': constitutive_call = 'yes'
                else: constitutive_call = 'no'
                ji.setConstitutive(constitutive_call)
                
                report_both_regions = 'no'
                try:
                    ### If the annotations are from a BED file produced by AltAnalyze, novel alternative splice sites may be present
                    ### if the below variable is not created, then this exon may over-ride the annotated exon region (e.g., E15.1 is over-written by E15.1_1234;E15.1_1256)
                    if 'ENS' in ji.JunctionID() and ':' not in ji.JunctionID(): report_both_regions = 'yes'
                except Exception: null=[]

                try:
                    ### If the annotations are from a BED file produced by AltAnalyze, it is possible for to a known exon to share a splice-site coordinate
                    ### with a novel junction exon. This will cause both to have the same override_annotation. Prevent this with the below 2nd override
                    if 'ENS' in ji.JunctionID() and ':' in ji.JunctionID(): override_annotation = None
                except Exception: null=[]
                
                if override_annotation != None:
                    if '_' in region1: region1 = string.split(region1,'_')[0]+'_'+str(int(string.split(region1,'_')[-1])-1)
                    if '_' in region2: region2 = string.split(region2,'_')[0]+'_'+str(int(string.split(region2,'_')[-1])+1)
                    if override_annotation == 1: region_id = region1 ### This forces a TopHat exon to be named for the splice-site position
                    else: region_id = region2                        
                else:
                    if report_both_regions == 'no':
                        ### Don't include specific start and end coordinates if inside a known exon
                        if ed1.AlignmentRegion() == 'exon': region1 = string.split(region1,'_')[0]
                        if ed2.AlignmentRegion() == 'exon': region2 = string.split(region2,'_')[0]
                        if ed1.AlignmentRegion() == 'full-intron' and ed2.AlignmentRegion() == 'full-intron':
                            region1 = string.split(region1,'_')[0]; region2 = string.split(region2,'_')[0]
                    
                    ### Below adjustmements need to compenstate for adjustments made upon import
                    if '_' in region1: region1 = string.split(region1,'_')[0]+'_'+str(int(string.split(region1,'_')[-1])-1)
                    if '_' in region2: region2 = string.split(region2,'_')[0]+'_'+str(int(string.split(region2,'_')[-1])+1)
                    
                ji.setExon1Stop(ji.Exon1Stop()-1); ji.setExon2Start(ji.Exon2Start()+1)
                if override_annotation != None: null=[] ### It is already assigned above
                elif region1 == region2: region_id = region1
                elif ji.Strand() == '+': region_id = region1+';'+region2
                else: region_id = region2+';'+region1 ### start and stop or genomically assigned
                uid = ji.GeneID()+':'+region_id
                #try: exon_region_db[ji.GeneID()].append((formatID(uid),region_id))
                #except KeyError: exon_region_db[ji.GeneID()]=[(formatID(uid),region_id)]
                ji.setExonRegionID(region_id)
                ji.setUniqueID(uid) ### hgu133
                ### Export format for new exons to add to the existing critical exon database (those in exon_region_db are combined with analyzed junctions)
                #exons_to_export[ji.GeneID(),region_id] = ji
            else:
                #print key, ji.GeneID(), ji.JunctionID(); sys.exit()
                null=[] ### Occurs because two genes are overlapping
    #return exons_to_export

def annotateNovelJunctions(novel_junction_db,novel_exon_db,exons_to_export):
    ### Reformat novel_exon_db to easily lookup exon data
    novel_exon_lookup_db={}
    for gene in novel_exon_db:
        for ed in novel_exon_db[gene]:
            novel_exon_lookup_db[gene,ed.ReadStart()]=ed
            
    ### Lookup the propper exon region ID and gene ID to format the unique ID and export coordinates
    junction_region_db={}
    for key in novel_junction_db:
        (chr,exon1_stop,exon2_start) = key
        ji=novel_junction_db[key]
        proceed = 'no'
        if ji.GeneID() != None:
            if ji.SpliceSitesFound() != 'both':
                e1 = (ji.GeneID(),exon1_stop)
                if ji.TransSplicing() == 'yes':
                    e2 = (ji.SecondaryGeneID(),exon2_start)
                else: e2 = (ji.GeneID(),exon2_start)
                if e1 in novel_exon_lookup_db and e2 in novel_exon_lookup_db:
                    proceed = 'yes'
                    try: ed1 = novel_exon_lookup_db[e1]; red1 = ed1.ExonRegionData(); gene1 = e1[0]
                    except Exception: print e1; kill
                    ed2 = novel_exon_lookup_db[e2]; red2 = ed2.ExonRegionData(); gene2 = e2[0]
                    ### If the splice-site was a match to a known junciton splice site, use it instead of that identified by exon-region location overlapp
                    if ji.LeftExonAnnotations() != None: region1 = ji.LeftExonAnnotations()
                    else: region1 = ed1.ExonRegionID(); exons_to_export[gene1,region1] = ed1
                    if ji.RightExonAnnotations() != None: region2 = ji.RightExonAnnotations()
                    else: region2 = ed2.ExonRegionID(); exons_to_export[gene2,region2] = ed2
                    #print region1,region2,ji.GeneID(),ji.Chr(),ji.Strand(), ji.LeftExonAnnotations(), ji.RightExonAnnotations()                    
            else:
                proceed = 'yes'
                region1 = ji.LeftExonAnnotations()
                region2 = ji.RightExonAnnotations()
                red1 = ji.LeftExonRegionData()
                red2 = ji.RightExonRegionData()
                ### Store the individual exons for export
                gene1 = ji.GeneID()
                if ji.TransSplicing() == 'yes': gene2 = ji.SecondaryGeneID()
                else: gene2 = ji.GeneID()
                exons_to_export[gene1,region1] = red1
                exons_to_export[gene2,region2] = red2
                
            if proceed == 'yes':
                try: splice_junctions = EnsemblImport.combineAnnotations([red1.AssociatedSplicingJunctions(),red2.AssociatedSplicingJunctions()])
                except Exception: print red1, red2;sys.exit()
                splice_events = EnsemblImport.combineAnnotations([red1.AssociatedSplicingEvent(),red2.AssociatedSplicingEvent()])
                ji.setAssociatedSplicingJunctions(splice_junctions)
                ji.setAssociatedSplicingEvent(splice_events)
                ens_exon_ids = EnsemblImport.combineAnnotations([red1.ExonID(),red2.ExonID()])
                ji.setExonID(ens_exon_ids)
                if ji.TransSplicing() == 'yes':
                    uid = ji.GeneID()+':'+region1+'-'+ji.SecondaryGeneID()+':'+region2
                    region_id = uid
                    ### When trans-splicing occurs, add the data twice to junction_region_db for the two different genes
                    ### in JunctionArray.inferJunctionComps, establish two separate gene junctions with a unique ID for the non-gene exon
                    try: junction_region_db[ji.GeneID()].append((formatID(uid),region1+'-'+'U1000.1_'+str(ji.Exon2Start())))
                    except KeyError: junction_region_db[ji.GeneID()]=[(formatID(uid),region1+'-'+'U1000.1_'+str(ji.Exon2Start()))]
                    try: junction_region_db[ji.SecondaryGeneID()].append((formatID(uid),'U0.1_'+str(ji.Exon1Stop())+'-'+region2))
                    except KeyError: junction_region_db[ji.SecondaryGeneID()]=[(formatID(uid),'U0.1_'+str(ji.Exon1Stop())+'-'+region2)]
                else:
                    uid = ji.GeneID()+':'+region1+'-'+region2
                    region_id = region1+'-'+region2
                    try: junction_region_db[ji.GeneID()].append((formatID(uid),region_id))
                    except KeyError: junction_region_db[ji.GeneID()]=[(formatID(uid),region_id)]
                ji.setExonRegionID(region_id)
                ji.setUniqueID(uid)
    return junction_region_db,exons_to_export

def alignReadsToExons(novel_exon_db,ens_exon_db,testImport=False):
    ### Simple method for aligning a single coordinate to an exon/intron region of an already matched gene
    examined_exons=0; aligned_exons=0
    for gene in ens_exon_db: #novel_exon_db
        try:
            region_numbers=[]; region_starts=[]; region_stops=[]
            for ed in novel_exon_db[gene]:
                examined_exons+=1; aligned_status=0; index=-1
                for rd in ens_exon_db[gene]:
                    index+=1 ### keep track of exon/intron we are in
                    region_numbers.append(int(string.split(rd.ExonRegionIDs()[1:],'.')[0]))
                    if rd.Strand() == '-': region_starts.append(rd.ExonStop()); region_stops.append(rd.ExonStart())
                    else: region_starts.append(rd.ExonStart()); region_stops.append(rd.ExonStop())
                    #print [rd.ExonStart(),rd.ExonStop(), rd.Strand()]
                    #print [ed.ReadStart(),rd.ExonStart(),rd.ExonStop()]
                    if ed.ReadStart()>=rd.ExonStart() and ed.ReadStart()<=rd.ExonStop():
                        ed.setAlignmentRegion('exon')
                        if 'I' in rd.ExonRegionIDs(): ### In an annotated intron
                            ed.setAlignmentRegion('intron')
                            ord = rd; updated = None
                            try: ### If the splice site is a novel 3' splice site then annotate as the 3' exon (less than 50nt away)
                                nrd = ens_exon_db[gene][index+1]
                                if (abs(ed.ReadStart()-nrd.ExonStart())<3) or (abs(ed.ReadStart()-nrd.ExonStop())<3):
                                    ed.setAlignmentRegion('full-intron') ### this is the start/end of intron coordinates
                                elif (abs(ed.ReadStart()-nrd.ExonStart())<50) or (abs(ed.ReadStart()-nrd.ExonStop())<50): rd = nrd; updated = 1
                            except Exception: null=[]
                            try:
                                prd = ens_exon_db[gene][index-1]
                                if (abs(ed.ReadStart()-prd.ExonStart())<3) or (abs(ed.ReadStart()-prd.ExonStop())<3):
                                    ed.setAlignmentRegion('full-intron')### this is the start/end of intron coordinates
                                elif (abs(ed.ReadStart()-prd.ExonStart())<50) or (abs(ed.ReadStart()-prd.ExonStop())<50):
                                    if updated==1: rd = ord; ###Hence the intron is too small to descriminate between alt5' and alt3' exons
                                    else: rd = prd
                            except Exception: null=[]                               
                        ed.setExonRegionData(rd); aligned_exons+=1; aligned_status=1
                        ed.setExonRegionID(rd.ExonRegionIDs()+'_'+str(ed.ReadStart()))
                        #print rd.ExonRegionIDs()+'_'+str(ed.ReadStart())
                        break
                if aligned_status == 0: ### non-exon/intron alinging sequences
                    region_numbers.sort(); region_starts.sort(); region_stops.sort()
                    if (rd.Strand() == '+' and ed.ReadStart()>=rd.ExonStop()) or (rd.Strand() == '-' and rd.ExonStop()>=ed.ReadStart()):
                        ### Applicable to 3'UTR (or other trans-splicing) aligning
                        utr_id = 'U'+str(region_numbers[-1])+'.1_'+str(ed.ReadStart())
                        ud = EnsemblImport.ExonAnnotationsSimple(rd.Chr(),rd.Strand(),region_stops[-1],region_stops[-1],gene,'','no',utr_id,'','')
                        ed.setExonRegionID(utr_id)
                    else:
                        ### Applicable to 5'UTR (or other trans-splicing) aligning
                        utr_id = 'U0.1'+'_'+str(ed.ReadStart())
                        ud = EnsemblImport.ExonAnnotationsSimple(rd.Chr(),rd.Strand(),region_starts[0],region_starts[0],gene,'','no',utr_id,'','')
                        ed.setExonRegionID(utr_id)
                    ed.setExonRegionData(ud)
                    ed.setAlignmentRegion('UTR')
        except Exception: null=[]
    if testImport == 'yes': print aligned_exons, 'splice sites aligned to exon region out of', examined_exons

def geneAlign(chr,chr_gene_locations,location_gene_db,chr_reads,switch_coord,read_aligned_to_gene):
    """ This function aligns the start or end position for each feature (junction or exon) to a gene, in two 
    steps by calling this function twice. In the second interation, the coordinates are reversed """
    
    index = 0 ### Don't examine genes already looked at
    genes_assigned = 0; trans_splicing=[]
    for (coord,ji) in chr_reads: ### junction coordinates or exon coordinates with gene object
        if index >5: index -=5 ### It is possible for some genes to overlap, so set back the index of genomically ranked genes each time
        gene_id_obtained = 'no'
        if switch_coord == 'no': rs,re=coord ### reverse the coordinates for the second iteration
        else: re,rs=coord ### first-interation coordinates (start and end)
        while index < len(chr_gene_locations):                 
            cs,ce = chr_gene_locations[index]
            ### Determine if the first listed coordinate lies within the gene
            if cs <= rs and ce >= rs:
                ### Yes, it does
                gene,strand = location_gene_db[chr,cs,ce]
                if switch_coord == 'yes': ### Only applies to coordinates, where the end-position didn't lie in the same gene as the start-position
                    if cs <= re and ce >= re:
                        ### This occurs when the first iteration detects a partial overlap, but the gene containing both coordinates is downstream
                        ### Hence, not trans-splicing
                        ji.setGeneID(gene)
                        break
                    first_geneid = ji.GeneID() ### see what gene was assigned in the first iteration (start position only)
                    #print ['trans',coord, first_geneid, gene]  ### Note: in rare cases, an exon can overlap with two genes (bad Ensembl annotations?)
                    ji.setTransSplicing()
                    side = ji.checkExonPosition(rs)
                    if side == 'left':
                        ji.setGeneID(gene)
                        ji.setSecondaryGeneID(first_geneid)
                    else:
                        ji.setSecondaryGeneID(gene)
                        #if ji.GeneID() == None: print 'B',coord, ji.GeneID(), secondaryGeneID()
                        #print ji.GeneID(), ji.SecondaryGeneID();kill
                    genes_assigned+=1; gene_id_obtained = 'yes'
                    ### Check to see if this gene represents a multi-gene spanning region (overlaps with multiple gene loci)
                    try:
                        cs2,ce2 = chr_gene_locations[index+1]
                        if cs2 < ce: index+=1 ### Continue analysis (if above is correct, the gene will have already been assigned)
                        else: break
                    except Exception: break
                else:
                    ### First iteration, store the identified gene ID (only looking at the start position)
                    ji.setGeneID(gene);  gene_id_obtained = 'yes'
                    #print gene, rs, re, cs, ce
                    ### Check the end position, to ensure it is also lies within the gene region
                    if cs <= re and ce >= re:
                        genes_assigned+=1
                    else:
                        ### Hence, the end lies outside the gene region
                        trans_splicing.append((coord,ji))
                    ### Check to see if this gene represents a multi-gene spanning region (overlaps with multiple gene loci)
                    try:
                        cs2,ce2 = chr_gene_locations[index+1]
                        if cs2 < ce: 
                            index+=1 ### Continue analysis (if above is correct, the gene will have already been assigned)
                        else: break
                    except Exception: break
            else:
                if rs < ce and re < ce: break
                elif switch_coord == 'no' and cs <= re and ce >= re:
                    ### This can occur if the left junction splice site is in an exon and the other is the UTR as opposed to another gene
                    gene,strand = location_gene_db[chr,cs,ce]
                    ji.setSecondaryGeneID(gene); gene_id_obtained = 'yes'
                    #print gene, coord, ji.Strand(), ji.GeneID()
                index+=1
        if gene_id_obtained == 'no':
            ### These often appear to be genes predicted by tBLASTn at UCSC but not by Ensembl (e.g., chr17:27,089,652-27,092,318 mouse mm9)
            null=[]
            #ji.setGeneID(None) ### This is not necessary, since if one exon does not align to a gene it is still a valid alignment
            #print chr,coord

    read_aligned_to_gene += genes_assigned           
    #print genes_assigned, chr, 'Gene IDs assigned out of', len(chr_reads)
    #print len(trans_splicing),'reads with evidence of trans-splicing'
    
    ### For any coordinate-pair where the end-position doesn't lie within the same gene as the start, re-run for those to see which gene they are in
    if switch_coord == 'no' and len(trans_splicing)>0:
        read_aligned_to_gene = geneAlign(chr,chr_gene_locations,location_gene_db,trans_splicing,'yes',read_aligned_to_gene)
    return read_aligned_to_gene


def getNovelExonCoordinates(species,root_dir):
    """ Currently, any novel exon determined during initial RNA-Seq read annotation with defined start and end coordinates, only has 
    the exon-end coordinate, not start, in it's name. However, the start and stop are indicated in the counts.Experiment.txt file.
    To get this, we parse that file and only store exons with an I or U in them and then correct for this in the matching function below """    
    exp_dir = root_dir+'/ExpressionInput/'
    dir_list = read_directory(exp_dir)
    counts_file = None
    for file in dir_list:
        if 'counts.' in file and 'steady' not in file:
            counts_file = file
    
    ### Example

    #ENSG00000137076:I17.1_35718353=chr9:35718353-35718403 (novel exon coordinates - just sorted, not necessarily in the correct order)
    #ENSG00000137076:E17.1-I17.1_35718403=chr9:35718809-35718403 (5' supporting junction)
    #ENSG00000137076:I17.1_35718353-E18.1=chr9:35718353-35717783 (3' supporting junction)
    #here, once we see that I17.1_35718353 is the exon ID, we know we need to get the function with -I17.1_35718403 (always the second value)
    
    if counts_file!=None:
        fn=filepath(exp_dir+counts_file)
        print 'Reading counts file'
        novel_exon_db = parseCountFile(fn,'exons',{}) ### Get novel exons
        print 'Reading counts file'
        novel_exon_db = parseCountFile(fn,'junctions',novel_exon_db) ### Get novel exons
    return novel_exon_db
    
def parseCountFile(fn,parseFeature,search_exon_db):
    novel_exon_db={}; firstLine=True
    for line in open(fn,'rU').xreadlines():
        key = string.split(line,'\t')[0]
        if firstLine: firstLine = False
        else:
            if '_' in key: ### Only look at novel exons
                #ENSG00000112695:I2.1_75953139=chr6:75953139-75953254
                uid, coordinates = string.split(key,'=')
                gene = string.split(uid,':')[0]
                if parseFeature == 'exons':
                    if '-' not in uid:
                        chr,coordinates = string.split(coordinates,':') ### Exclude the chromosome
                        coord1,coord2 = string.split(coordinates,'-')
                        intron = string.split(uid,'_')[0]
                        intron = string.split(intron,':')[1]
                        first = intron+'_'+coord1
                        second = intron+'_'+coord2
                        proceed = True
                        if first in uid: search_uid = second ### if the first ID is already the one looked for, store the second with the exon ID
                        elif second in uid: search_uid = first
                        else:
                            proceed = False
                            #print uid, first, second; sys.exit()
                            #example: ENSG00000160785:E2.15_156170151;E2.16_156170178=chr1:156170151-156170178
                        if proceed:
                            try: novel_exon_db[gene].append((uid,search_uid))
                            except Exception: novel_exon_db[gene] = [(uid,search_uid)]  
                elif '-' in uid and 'I' in uid: ### get junctions
                    if gene in search_exon_db:
                        for (u,search_uid) in search_exon_db[gene]:
                            #if gene == 'ENSG00000137076': print u,search_uid,uid
                            if search_uid in uid:
                                novel_exon_db[uid] = u ### Relate the currently examined novel exon ID to the junction not current associated
                                #if gene == 'ENSG00000137076': print u, uid
                                #print uid;sys.exit()
    return novel_exon_db

def getHighExpNovelExons(fn):
    """ Idea - if the ranking of exons based on expression changes from one condition to another, alternative splicing is occuring """
    exon_max_exp_db={}; uid_key_db={}; firstLine=True
    cutoff = 0.5
    for line in open(fn,'rU').xreadlines():
        t = string.split(line,'\t')
        if firstLine: firstLine = False
        else:
            key=t[0]
            #ENSG00000112695:I2.1_75953139=chr6:75953139-75953254
            uid, coordinates = string.split(key,'=')
            gene = string.split(uid,':')[0]
            max_read_counts = max(map(lambda x: float(x), t[1:]))
            try: exon_max_exp_db[gene].append((max_read_counts,uid))
            except Exception: exon_max_exp_db[gene] = [(max_read_counts,uid)]
            uid_key_db[uid] = key ### retain the coordinate info
    
    out_file = string.replace(fn,'.txt','-highExp.txt')
    out_obj = export.ExportFile(out_file)
    ### Compare the relative expression of junctions and exons separately for each gene (junctions are more comparable)
    for gene in exon_max_exp_db:
        junction_set=[]; exon_set=[]
        exon_max_exp_db[gene].sort()
        exon_max_exp_db[gene].reverse()
        for (exp,uid) in exon_max_exp_db[gene]:
            if '-' in uid: junction_set.append((exp,uid))
            else: exon_set.append((exp,uid))
        if len(junction_set)>0:
            maxJunctionExp = junction_set[0][0]
            junction_percent_exp = map(lambda x: (x[1],expThreshold(x[0]/maxJunctionExp,cutoff)), junction_set)
            high_exp_junctions = []
            for (uid,p) in junction_percent_exp: ### ID and percentage of expression
                if p!='NA':
                    out_obj.write(uid_key_db[uid]+'\t'+p+'\n') ### write out the original ID with coordinates
            
        if len(exon_set)>0:
            maxExonExp = exon_set[0][0]
            exon_percent_exp = map(lambda x: (x[1],expThreshold(x[0]/maxExonExp,cutoff)), exon_set)
            high_exp_exons = []
            for (uid,p) in exon_percent_exp: ### ID and percentage of expression
                if p!='NA':
                    out_obj.write(uid_key_db[uid]+'\t'+p+'\n') 
    out_obj.close()


def expThreshold(ratio,cutoff):
    #print [ratio,cutoff]
    if ratio>cutoff: return str(ratio)
    else: return 'NA'

def compareExonAndJunctionResults(species,array_type,summary_results_db,root_dir):
    results_dir = root_dir +'AltResults/AlternativeOutput/'
    dir_list = read_directory(results_dir)
    filtered_dir_db={}
    try: novel_exon_junction_db = getNovelExonCoordinates(species,root_dir)
    except Exception:
        print 'No counts file found.'
        novel_exon_junction_db={} ### only relevant to RNA-Seq analyses
    for comparison_file in summary_results_db:
        for results_file in dir_list:
            if (comparison_file in results_file and '-exon-inclusion-results.txt' in results_file) and ('comparison' not in results_file):
                try: filtered_dir_db[comparison_file].append(results_file)
                except Exception: filtered_dir_db[comparison_file] = [results_file]
    for comparison_file in filtered_dir_db:
        alt_result_files = filtered_dir_db[comparison_file]
        #print alt_result_files, comparison_file
        importAltAnalyzeExonResults(alt_result_files,novel_exon_junction_db,results_dir)
        
    ### Build combined clusters of high-confidence exons
    import ExpressionBuilder
    try:
        input_dir = string.split(results_dir,'AltResults')[0]+'GO-Elite/AltExonConfirmed/'
        cluster_file, rows_in_file = ExpressionBuilder.buildAltExonClusterInputs(input_dir,species,array_type,dataType='AltExonConfirmed')
        if rows_in_file > 8000: useHOPACH = False
        else: useHOPACH = True
        if rows_in_file < 9000:
            ExpressionBuilder.exportHeatmap(cluster_file)
    except Exception: pass

    try:
        input_dir = string.split(results_dir,'AltResults')[0]+'GO-Elite/AltExon/'
        cluster_file, rows_in_file = ExpressionBuilder.buildAltExonClusterInputs(input_dir,species,array_type,dataType='AltExon')
        if rows_in_file > 8000: useHOPACH = False
        else: useHOPACH = True
        if rows_in_file < 9000:
            ExpressionBuilder.exportHeatmap(cluster_file)
    except Exception: pass

class SplicingData:
    def __init__(self,score,symbol,description,exonid,probesets,direction,splicing_event,external_exon,genomic_loc,gene_exp,protein_annot,domain_inferred,domain_overlap,method,dataset):
        self.score = score; self.dataset = dataset
        self.symbol = symbol;
        self.description=description;self.exonid=exonid;self.probesets=probesets;self.direction=direction
        self.splicing_event=splicing_event;self.external_exon=external_exon;self.genomic_loc=genomic_loc;
        self.gene_exp=gene_exp;self.protein_annot=protein_annot;self.domain_inferred=domain_inferred
        self.domain_overlap=domain_overlap;self.method=method
    def Score(self): return self.score
    def setScore(self,score): self.score = score
    def GeneExpression(self): return self.gene_exp
    def Dataset(self): return self.dataset
    def Symbol(self): return self.symbol
    def Description(self): return self.description
    def ExonID(self): return self.exonid
    def appendExonID(self,exonid): self.exonid+='|'+exonid
    def Probesets(self): return self.probesets
    def ProbesetDisplay(self):
        if len(self.Probesets()[1])>0:
            return string.join(self.Probesets(),'-')
        else:
            return self.Probesets()[0]
    def ProbesetsSorted(self):
        ### Don't sort the original list
        a = [self.probesets[0],self.probesets[1]]
        a.sort()
        return a
    def Direction(self): return self.direction
    def setDirection(self,direction): self.direction = direction
    def SplicingEvent(self): return self.splicing_event
    def ProteinAnnotation(self): return self.protein_annot
    def DomainInferred(self): return self.domain_inferred
    def DomainOverlap(self): return self.domain_overlap
    def Method(self): return self.method
    def setEvidence(self,evidence): self.evidence = evidence
    def Evidence(self): return self.evidence
    def GenomicLocation(self): return self.genomic_loc
    def setExonExpStatus(self, exon_expressed): self.exon_expressed = exon_expressed
    def ExonExpStatus(self): return self.exon_expressed
    
def importAltAnalyzeExonResults(dir_list,novel_exon_junction_db,results_dir):
    
    regulated_critical_exons={}; converted_db={}
    includeExonJunctionComps=True ### Allow ASPIRE comparisons with the inclusion feature as an exon to count for additive reciprocal evidence
    print "Reading AltAnalyze results file"
    root_dir = string.split(results_dir,'AltResults')[0]
    
    for filename in dir_list:
        x=0; regulated_critical_exon_temp={}
        fn=filepath(results_dir+filename)
        new_filename = string.join(string.split(filename,'-')[:-5],'-')
        if '_vs_' in filename and '_vs_' in new_filename: export_filename = new_filename
        else: export_filename = string.join(string.split(filename,'-')[:-5],'-')

        export_path = results_dir+export_filename+'-comparison-evidence.txt'
        try: os.remove(filepath(export_path)) ### If we don't do this, the old results get added to the new
        except Exception: null=[]
    
        if 'AltMouse' in filename:
            altmouse_ensembl_db = importAltMouseEnsembl()
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if x==0: x=1; #print t[12],t[13],t[22],t[23]
            else:
                converted = False ### Indicates both junction sides were regulated
                geneid = t[0]; exonid = t[4]; probeset1 = t[6]; probeset2 = ''; score = t[1][:4]; symbol = t[2]; description = t[3]; regions = t[-4]; direction = t[5]
                genomic_loc = t[-1]; splicing_event = t[-3]; external_exon = t[-6]; gene_exp_fold = t[-8]; protein_annot = t[14]; domain_inferred = t[15]; domain_overlap = t[17]
                expressed_exon = 'NA'
                if 'RNASeq' in filename: expressed_exon = 'no' ### Set by default
                if ':' in geneid: geneid = string.split(geneid,':')[0] ### User reported that gene:gene was appearing and not sure exactly where or why but added this to address it
                if 'FIRMA' in fn: method = 'FIRMA'
                elif 'splicing-index' in fn: method = 'splicing-index'
                if 'ASPIRE' in filename or 'linearregres' in filename:
                    f1=float(t[12]); f2=float(t[13]); probeset1 = t[8]; probeset2 = t[10]; direction = t[6]; exonid2 = t[5]; splicing_event = t[-4]
                    protein_annot = t[19]; domain_inferred = t[20]; domain_overlap = t[24]; method = 'linearregres'; regions = t[-5]
                    exon1_exp=float(t[-15]); exon2_exp=float(t[-14]); fold1=float(t[12]); fold2=float(t[13])
                    if fold1<0: fold1 = 1 ### don't factor in negative changes
                    if fold2<0: fold2 = 1 ### don't factor in negative changes
                    """
                    if 'RNASeq' not in filename:
                        exon1_exp = math.pow(2,exon1_exp)
                        exon2_exp = math.log(2,exon2_exp)
                    
                    m1 = exon1_exp*fold1
                    m2 = exon2_exp*fold2
                    max_exp = max([m1,m2])
                    min_exp = min([m1,m2])
                    percent_exon_expression = str(min_exp/max_exp)
                    """
                    if 'ASPIRE' in filename: method = 'ASPIRE'; score = t[1][:5]
                    if '-' not in exonid and includeExonJunctionComps == False:
                        exonid=None ### Occurs when the inclusion just in an exon (possibly won't indicate confirmation so exclude)
                    else: exonid = exonid+' vs. '+exonid2
                    if 'AltMouse' in filename:
                        try: geneid = altmouse_ensembl_db[geneid]
                        except Exception: geneid = geneid
                    if 'RNASeq' not in filename and 'junction' not in filename: regions = string.replace(regions,'-','.')
                else:
                    if 'RNASeq' in filename and '-' not in exonid:
                        fold = float(t[10]); exon_exp = float(t[18]); gene_exp = float(t[19])
                        if fold < 0: fold = -1.0/fold
                        GE_fold = float(gene_exp_fold)
                        if GE_fold < 0: GE_fold = -1.0/float(gene_exp_fold)
                        exon_psi1 = abs(exon_exp)/(abs(gene_exp))
                        exon_psi2 = (abs(exon_exp)*fold)/(abs(gene_exp)*GE_fold)
                        max_incl_exon_exp = max([exon_psi1,exon_psi2])
                        #if max_incl_exon_exp>0.20: expressed_exon = 'yes'
                        expressed_exon = max_incl_exon_exp
                        #if 'I2.1_75953139' in probeset1:
                        #print [exon_exp,gene_exp,exon_exp*fold,gene_exp*GE_fold]
                        #print exon_psi1, exon_psi2;sys.exit()
                probesets = [probeset1,probeset2]
                if (method == 'splicing-index' or method == 'FIRMA') and ('-' in exonid) or exonid == None:
                    pass #exclude junction IDs
                else:
                    regions = string.replace(regions,';','|')
                    regions = string.replace(regions,'-','|')
                    regions = string.split(regions,'|')
                    for region in regions:
                        if len(region) == 0:
                            try: region = t[17]+t[18] ### For junction introns where no region ID exists
                            except Exception: null=[]
                        if ':' in region: region = string.split(region,':')[-1] ### User reported that gene:gene was appearing and not sure exactly where or why but added this to address it
                        if probeset1 in novel_exon_junction_db:
                            uid = novel_exon_junction_db[probeset1] ### convert the uid (alternative exon) to the annotated ID for the novel exon
                            converted_db[uid] = probeset1
                        else:
                            uid = geneid+':'+region
                        ss = SplicingData(score,symbol,description,exonid,probesets,direction,splicing_event,external_exon,genomic_loc,gene_exp_fold,protein_annot,domain_inferred,domain_overlap,method,filename)
                        ss.setExonExpStatus(str(expressed_exon))
                        try: regulated_critical_exon_temp[uid].append(ss)
                        except Exception: regulated_critical_exon_temp[uid] = [ss]

        
        #print filename, len(regulated_critical_exon_temp)
        for uid in regulated_critical_exon_temp:
            report=None
            if len(regulated_critical_exon_temp[uid])>1:
                ### We are only reporting one here and that's OK, since we are only reporting the top scores... won't include all inclusion junctions.
                scores=[]
                for ss in regulated_critical_exon_temp[uid]: scores.append((float(ss.Score()),ss))
                scores.sort()
                if (scores[0][0]*scores[-1][0])<0:
                    ss1 = scores[0][1]; ss2 = scores[-1][1]
                    if ss1.ProbesetsSorted() == ss2.ProbesetsSorted(): ss1.setDirection('mutual') ### same exons, hence, mutually exclusive event (or similiar)
                    else: ss1.setDirection('both') ### opposite directions in the same comparison-file, hence, conflicting data
                    report=[ss1]
                else:
                    if abs(scores[0][0])>abs(scores[-1][0]): report=[scores[0][1]]
                    else: report=[scores[-1][1]]
            else:
                report=regulated_critical_exon_temp[uid]
            ### Combine data from different analysis files
            try: regulated_critical_exons[uid]+=report
            except Exception: regulated_critical_exons[uid]=report
            """if 'ENSG00000204120' in uid:
                print uid,
                for i in regulated_critical_exon_temp[uid]:
                    print i.Probesets(),
                print ''
                """
            try: report[0].setEvidence(len(regulated_critical_exon_temp[uid])) ###set the number of exons demonstrating regulation of this exons
            except Exception: null=[]

    clearObjectsFromMemory(regulated_critical_exon_temp)
    
    export_data,status = AppendOrWrite(export_path)
    
    if status == 'not found': 
        header = string.join(['uid','source-IDs','symbol','description','exonids','independent confirmation','score','regulation direction','alternative exon annotations','associated isoforms','inferred regulated domains','overlapping domains','method','supporting evidence score','novel exon: high-confidence','percent exon expression of gene','differential gene-expression','genomic location'],'\t')+'\n'
        export_data.write(header)

    print len(regulated_critical_exons), 'regulated exon IDs imported.\n'
    print 'writing:',export_path; n=0
    
    ### Check for alternative 3' or alternative 5' exon regions that were not matched to the right reciprocal junctions (occurs because only one of the exon regions is called alternative)
    regulated_critical_exons_copy={}
    for uid in regulated_critical_exons:
        regulated_critical_exons_copy[uid]=regulated_critical_exons[uid]
    u=0
    ### This is most applicable to RNA-Seq since the junction IDs correspond to the Exon Regions not the probeset Exon IDs
    for uid in regulated_critical_exons_copy: ### Look through the copied version since we can't delete entries while iterating through
        ls = regulated_critical_exons_copy[uid]
        u+=1
        #if u<20: print uid
        for jd in ls:
            if jd.Method() != 'splicing-index' and jd.Method() != 'FIRMA':
                try: ### Applicable to RNA-Seq
                    gene,exonsEx = string.split(jd.Probesets()[1],':') ### Exclusion probeset will have the exon not annotated as the critical exon (although it should be as well)
                    gene,exonsIn = string.split(jd.Probesets()[0],':')
                except Exception:
                    gene, ce = string.split(uid,':')
                    exonsIn, exonsEx = string.split(jd.ExonID(),'vs.')
                if gene !=None:
                    critical_exon = None
                    five_prime,three_prime = string.split(exonsEx,'-')
                    try: five_primeIn,three_primeIn = string.split(exonsIn,'-')
                    except Exception: five_primeIn = exonsIn; three_primeIn = exonsIn ### Only should occur during testing when a exon rather than junction ID is considered
                    #if gene == 'ENSG00000133083': print five_prime,three_prime, five_primeIn,three_primeIn
                    if five_primeIn == five_prime: ### Hence, the exclusion 3' exon should be added
                        critical_exon = gene+':'+three_prime
                        exonid = three_prime
                    elif three_primeIn == three_prime: ### Hence, the exclusion 3' exon should be added
                        critical_exon = gene+':'+five_prime
                        exonid = five_prime
                    else:
                        if ('5' in jd.SplicingEvent()) or ('five' in jd.SplicingEvent()):
                            critical_exon = gene+':'+five_prime
                            exonid = five_prime
                        elif ('3' in jd.SplicingEvent()) or ('three' in jd.SplicingEvent()):
                            critical_exon = gene+':'+three_prime
                            exonid = three_prime
                        elif ('alt-N-term' in jd.SplicingEvent()) or ('altPromoter' in jd.SplicingEvent()):
                            critical_exon = gene+':'+five_prime
                            exonid = five_prime
                        elif ('alt-C-term' in jd.SplicingEvent()):
                            critical_exon = gene+':'+three_prime
                            exonid = three_prime
                    #print critical_exon, uid, jd.ExonID(),jd.SplicingEvent(); sys.exit() 
                    if critical_exon != None:
                        if critical_exon in regulated_critical_exons:
                            #print uid, critical_exon; sys.exit()
                            if len(regulated_critical_exons[critical_exon]) == 1:
                                if len(ls)==1 and uid in regulated_critical_exons: ### Can be deleted by this method
                                    if 'vs.' not in regulated_critical_exons[critical_exon][0].ExonID() and 'vs.' not in regulated_critical_exons[critical_exon][0].ExonID():
                                        regulated_critical_exons[uid].append(regulated_critical_exons[critical_exon][0])
                                        del regulated_critical_exons[critical_exon]
                                elif uid in regulated_critical_exons: ###If two entries already exit
                                    ed = regulated_critical_exons[uid][1]
                                    ed2 = regulated_critical_exons[critical_exon][0]
                                    if 'vs.' not in ed.ExonID() and 'vs.' not in ed2.ExonID():
                                        if ed.Direction() != ed2.Direction(): ### should be opposite directions
                                            ed.appendExonID(exonid)
                                            ed.setEvidence(ed.Evidence()+1)
                                            ed.setScore(ed.Score()+'|'+ed2.Score())
                                            del regulated_critical_exons[critical_exon]
                                            
    firstEntry=True
    for uid in regulated_critical_exons:
        if uid in converted_db:
            converted = True
        else: converted = False
        #if 'ENSG00000133083' in uid: print [uid]
        exon_level_confirmation = 'no'
        ls = regulated_critical_exons[uid]
        jd = regulated_critical_exons[uid][0] ### We are only reporting one here and that's OK, since we are only reporting the top scores... won't include all inclusion junctions.
        if len(ls)>1:
            methods = []; scores = []; direction = []; exonids = []; probesets = []; evidence = 0; genomic_location = []
            junctionids=[]
            junction_data_found = 'no'; exon_data_found = 'no'
            for jd in ls:
                if jd.Method() == 'ASPIRE' or jd.Method() == 'linearregres': 
                    junction_data_found = 'yes'
                    methods.append(jd.Method())
                    scores.append(jd.Score())
                    direction.append(jd.Direction())
                    exonids.append(jd.ExonID())
                    junctionids.append(jd.ExonID())
                    probesets.append(jd.ProbesetDisplay())
                    evidence+=jd.Evidence()
                    genomic_location.append(jd.GenomicLocation())
                    ### Prefferentially obtain isoform annotations from the reciprocal analysis which is likely more accurate
                    isoform_annotations = [jd.ProteinAnnotation(), jd.DomainInferred(), jd.DomainOverlap()]
            for ed in ls:
                if ed.Method() == 'splicing-index' or ed.Method() == 'FIRMA':
                    exon_data_found = 'yes' ### pick one of them
                    methods.append(ed.Method())
                    scores.append(ed.Score())
                    direction.append(ed.Direction())
                    exonids.append(ed.ExonID())
                    probesets.append(ed.ProbesetDisplay())
                    evidence+=ed.Evidence()
                    genomic_location.append(ed.GenomicLocation())
                    #isoform_annotations = [ed.ProteinAnnotation(), ed.DomainInferred(), ed.DomainOverlap()]
            if junction_data_found == 'yes' and exon_data_found == 'yes':
                exon_level_confirmation = 'yes'
                for junctions in junctionids:
                    if 'vs.' in junctions:
                        j1 = string.split(junctions,' vs. ')[0] ### inclusion exon or junction
                        if '-' not in j1: ### not a junction, hence, may not be sufficient to use for confirmation (see below)
                            if 'I' in j1: ### intron feature
                                if '_' in j1: ### novel predicted exon
                                    exon_level_confirmation = 'no'
                                else:
                                    exon_level_confirmation = 'yes'
                            else:
                                if '_' in j1:
                                    exon_level_confirmation = 'no'
                                else:
                                    exon_level_confirmation = 'partial'
            method = string.join(methods,'|')
            unique_direction = unique.unique(direction)
            genomic_location = unique.unique(genomic_location)
            if len(unique_direction) == 1: direction = unique_direction[0]
            else: direction = string.join(direction,'|')
            score = string.join(scores,'|')
            probesets = string.join(probesets,'|')
            exonids_unique = unique.unique(exonids)
            if len(exonids_unique) == 1: exonids = exonids_unique[0]
            else: exonids = string.join(exonids,'|')
            if len(genomic_location) == 1: genomic_location = genomic_location[0]
            else: genomic_location = string.join(genomic_location,'|')
            evidence = str(evidence)
            if 'mutual' in direction: direction = 'mutual'
        if len(ls) == 1:
            probesets = jd.ProbesetDisplay()
            direction = jd.Direction()
            score = jd.Score()
            method = jd.Method()
            exonids = jd.ExonID()
            evidence = jd.Evidence()
            genomic_location = jd.GenomicLocation()
            isoform_annotations = [jd.ProteinAnnotation(), jd.DomainInferred(), jd.DomainOverlap()]
        try:
            if int(evidence)>4 and 'I' in uid: novel_exon = 'yes' ### high-evidence novel exon
            else: novel_exon = 'no'
            if converted == True:novel_exon = 'yes'
            else: novel_exon = 'no'
            values = [uid, probesets, jd.Symbol(), jd.Description(), exonids, exon_level_confirmation, score, direction, jd.SplicingEvent()]
            values += isoform_annotations+[method, str(evidence),novel_exon,jd.ExonExpStatus(),jd.GeneExpression(),genomic_location]
            values = string.join(values,'\t')+'\n'
            #if 'yes' in exon_level_confirmation:
            export_data.write(values); n+=1
            if exon_level_confirmation == 'yes' and ('|' not in direction):
                geneID = string.split(uid,':')[0]
                if firstEntry:
                    ### Also export high-confidence predictions for GO-Elite
                    elite_export_path = string.split(results_dir,'AltResults')[0]+'GO-Elite/AltExonConfirmed/'+export_filename+'-junction-exon-evidence2.txt'
                    elite_export_path = string.split(results_dir,'AltResults')[0]+'GO-Elite/AltExonConfirmed/junction-exon-evidence2.txt'
                    #elite_export_data = export.ExportFile(elite_export_path)
                    #elite_export_data.write('GeneID\tEn\tExonID\tScores\tGenomicLocation\n')
                    elite_export_data = AppendOrWrite(elite_export_path)[0]
                    firstEntry = False
                #if 'DNA' in isoform_annotations[-1]:
                if '2moter' not in jd.SplicingEvent() and '2lt-N' not in jd.SplicingEvent():
                        values = [uid, probesets, jd.Symbol(), jd.Description(), exonids, exon_level_confirmation, score, direction, jd.SplicingEvent()]
                        values += isoform_annotations+[method, str(evidence),novel_exon,jd.ExonExpStatus(),jd.GeneExpression(),genomic_location,export_filename]
                        values = string.join(values,'\t')+'\n'
                        elite_export_data.write(values)

        except Exception:
            pass ### Unknown error - not evaluated in 2.0.8  - isoform_annotations not referenced
    print n,'exon IDs written to file.'
    export_data.close()
    try: elite_export_data.close()
    except Exception: pass
    clearObjectsFromMemory(regulated_critical_exons)
    clearObjectsFromMemory(regulated_critical_exons_copy)
    #print '!!!!Within comparison evidence'
    #returnLargeGlobalVars()

if __name__ == '__main__':
    #reformatExonFile('Hs','exon',True);sys.exit()
    filename = '/Volumes/Time Machine Backups/dataAnalysis/PCBC_Sep2013/C4-reference/ExpressionInput/counts.C4.txt'
    #fastRPKMCalculate(filename);sys.exit()
    file1 = '/Volumes/Time Machine Backups/dataAnalysis/CardiacRNASeq/versus D0/ExpressionInput/counts.CardiacRNASeq.txt'
    file2 = '/Volumes/Time Machine Backups/dataAnalysis/PCBC_Sep2013/C4-bedFiles/ExpressionInput/counts.C4.txt'
    #getHighExpNovelExons(file1);sys.exit()
    #mergeCountFiles(file1,file2); sys.exit()
    import UI
    test_status = 'yes'
    data_type = 'ncRNA'
    data_type = 'mRNA'
    array_type = 'RNASeq'
    species = 'Hs' ### edit this
    
    summary_results_db = {}
    root_dir = '/Volumes/Time Machine Backups 1/dataAnalysis/CardiacRNASeq/StringentAnalyses/'
    root_dir = '/Volumes/My Passport/dataAnalysis/CardiacRNASeq/BedFiles/'
    #root_dir = '/Volumes/Time Machine Backups/dataAnalysis/PCBC_Sep2013/C4-COI_treatment/'
    #summary_results_db['Hs_Junction_d14_vs_d7.p5_average-ASPIRE-exon-inclusion-results.txt'] = [] ### edit this
    #summary_results_db['Hs_Junction_d14_vs_d7.p5_average-splicing-index-exon-inclusion-results.txt'] = [] ### edit this
    
    results_dir = root_dir +'AltResults/AlternativeOutput/'
    dir_list = read_directory(results_dir)
    for i in dir_list:
        if '_average' in i:
            comparison, end = string.split(i,'_average')
            if '-exon-inclusion-results.txt' in i: summary_results_db[comparison]=[]

    compareExonAndJunctionResults(species,array_type,summary_results_db,root_dir); sys.exit()
    
    fl = UI.ExpressionFileLocationData('','','',''); fl.setCELFileDir(loc); fl.setRootDir(loc)
    exp_file_location_db={}; exp_file_location_db['test']=fl
    alignJunctionsToEnsembl(species,exp_file_location_db,'test'); sys.exit()
    getEnsemblAssociations(species,data_type,test_status,'yes'); sys.exit()
