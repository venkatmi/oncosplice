#!/usr/bin/env python

#!/usr/bin/env python
import numpy as np
import sys,string
import os
import os.path
from collections import defaultdict
import sampleIndexSelection
import export
from bs4 import BeautifulSoup
import urllib2
import csv
import sys
import re
import math
import export


def cell_text(cell):
    return " ".join(cell.stripped_strings)


    
def parseResultfolders(motifdir,GEdir,SFlist):
    sfs=[]
    for lin in open(SFlist,'rU').xreadlines():
        s=lin.rstrip('\r\n')
        s1=string.split(s,'\t')
        sfs.append(s1[0])
    
    mappingdict=defaultdict(list)
    allden=[]
    for filename in os.listdir(motifdir):
     name=filename
     mapping=[]
     dellst=[]
     
     if "._" not in filename and "Events" not in filename:
        fol=os.path.join(motifdir,filename)
        if os.path.isdir(fol):
            #for filename2 in os.listdir(fol):
             #filnam2=os.path.join(fol,filename2)
             #if "._" not in filnam2:
             #   if os.path.isdir(filnam2):
             #       #print filnam2
             #       flag=0
             #       if "._" not in filename2:
             #           name=filename+":"+filename2
             #           flag=1
             #       
             #       if flag==1:
                        for filename3 in os.listdir(fol):
                           
                            
                            if filename3=="finalResults.tab":
                                
                                clipres=os.path.join(fol,filename3)
                                for lin in open(clipres,'rU').xreadlines():
                                
                                    q=lin.rstrip('\r\n')
                                    q1=string.split(q,'\t')

                                    
                                    clipnam=q1[0]+":"+q1[1]+":"+q1[2]
                                    mappingdict[name,clipnam,"Clipseq"]=q1[11]
                           
                            if filename3=="output_TF_strand":
                                knownrbp=os.path.join(fol,filename3)
                                for filename4 in os.listdir(knownrbp):
                                    if filename4=="knownResults.txt":
                                        filenam4=os.path.join(knownrbp,filename4)
                                        try:
                                            head=0
                                            for line in open(filenam4,'rU').xreadlines():
                                                q=line.rstrip('\r\n')
                                                q1=string.split(q,'\t')
                                                if head==0:
                                                    motif=q1.index('Motif Name')
                                                    pval=q1.index('P-value')
                                                    head=1
                                                    continue
                                                else:
                                                    mappingdict[name,q1[motif],"Cisbp_Actual"]=q1[pval]
                                                    
                                        except Exception:
                                            continue
                            
                            if filename3=="output1":
                                denovorbp=os.path.join(fol,filename3)
                                for filename4 in os.listdir(denovorbp):
                                    if filename4=="homerResults.html":
                                        denolink="file://"+str(os.path.join(denovorbp,filename4))
                                        #print denolink
                                        html = urllib2.urlopen(denolink).read()
                                        soup = BeautifulSoup(html)
                                        for table in soup.find_all('table'):
                                            for row in table.find_all('tr'):
                                                col = map(cell_text, row.find_all(re.compile('t[dh]')))
                                                
                                                
                                                if col[2]=="P-value":
                                                    continue
                                                else:

                                                    motname=string.split(col[7],"(")[0]
                                                    mapping.append([name+";"+motname,float(col[2])])
                                                    #mappingdict[name,motname,"Cisbp_denovo"]=col[2]
                            
                            if filename3=="output2":
                                denovorbp=os.path.join(fol,filename3)
                                for filename4 in os.listdir(denovorbp):
                                    if filename4=="homerResults.html":
                                        denolink="file://"+str(os.path.join(denovorbp,filename4))
                                        #print denolink
                                        html = urllib2.urlopen(denolink).read()
                                        soup = BeautifulSoup(html)
                                        for table in soup.find_all('table'):
                                            for row in table.find_all('tr'):
                                                col = map(cell_text, row.find_all(re.compile('t[dh]')))
                                                if col[2]=="P-value":
                                                    continue
                                                else:
                                                    motname=string.split(col[7],"(")[0]
                                                    mapping.append([name+";"+motname,float(col[2])])
                                                    #mappingdict[name,motname,"Cisbp_denovo"]=col[2]
                            if filename3=="output3":
                                denovorbp=os.path.join(fol,filename3)
                                for filename4 in os.listdir(denovorbp):
                                    if filename4=="homerResults.html":
                                        denolink="file://"+str(os.path.join(denovorbp,filename4))
                                        #print denolink
                                        html = urllib2.urlopen(denolink).read()
                                        soup = BeautifulSoup(html)
                                        for table in soup.find_all('table'):
                                            for row in table.find_all('tr'):
                                                col = map(cell_text, row.find_all(re.compile('t[dh]')))
                                                if col[2]=="P-value":
                                                    continue
                                                else:
                                                    motname=string.split(col[7],"(")[0]
                                                    mapping.append([name+";"+motname,float(col[2])])
                                                    #mappingdict[name,motname,"Cisbp_denovo"]=col[2]
                            if filename3=="output4":
                                denovorbp=os.path.join(fol,filename3)
                                for filename4 in os.listdir(denovorbp):
                                    if filename4=="homerResults.html":
                                        denolink="file://"+str(os.path.join(denovorbp,filename4))
                                        #print denolink
                                        html = urllib2.urlopen(denolink).read()
                                        soup = BeautifulSoup(html)
                                        for table in soup.find_all('table'):
                                            for row in table.find_all('tr'):
                                                col = map(cell_text, row.find_all(re.compile('t[dh]')))
                                                if col[2]=="P-value":
                                                    continue
                                                else:
                                                    motname=string.split(col[7],"(")[0]
                                                    mapping.append([name+";"+motname,float(col[2])])
                                                    #mappingdict[name,motname,"Cisbp_denovo"]=col[2]
                            if filename3=="output5":
                                denovorbp=os.path.join(fol,filename3)
                                for filename4 in os.listdir(denovorbp):
                                    if filename4=="homerResults.html":
                                        denolink="file://"+str(os.path.join(denovorbp,filename4))
                                        #print denolink
                                        html = urllib2.urlopen(denolink).read()
                                        soup = BeautifulSoup(html)
                                        for table in soup.find_all('table'):
                                            for row in table.find_all('tr'):
                                                col = map(cell_text, row.find_all(re.compile('t[dh]')))
                                                if col[2]=="P-value":
                                                    continue
                                                else:
                                                    motname=string.split(col[7],"(")[0]
                                                    mapping.append([name+";"+motname,float(col[2])])
                                                    #print name,motname,col[2]
                                                    #sys.exit()
                                                    #mappingdict[name,motname,"Cisbp_denovo"]=col[2]
     mapping.sort(key=lambda x: x[0])
                            
     mapping.sort(key=lambda x: x[1])
                            #prev=""
     #output=os.path.join(motifdir,"test.txt")
     #output_w=open(output,"a")
     for i in range(len(mapping)):
        if mapping[i][0] not in dellst:
            mot=string.split(mapping[i][0],";")[1]
            genes=[]
            genes=string.split(mot,":")[1:]
            allden.append([filename,mot,genes,mapping[i][1]])
            #output_w.write(mapping[i][0]+"\t"+str(mapping[i][1]))
      #      output_w.write("\n")
            dellst.append(mapping[i][0])
     final={}                       
     for i in range(len(allden)):
        de=[]
        de= allden[i]
        
        for q in de[2]:
                if q in final:
                    if de[3] < final[q][1]:
                        final[q]=[de[0],de[3],de[1]]
                else:
                    final[q]=[de[0],de[3],de[1]]
     for genes in final:
      
      de=[]
      de=final[genes]
      mappingdict[de[0],de[2],"Cisbp_denovo"]=str(de[1])
     

    for filename in os.listdir(GEdir):
        if "GE" in filename and "._GE" not in filename:
            InputFile=os.path.join(GEdir, filename)
            name=string.replace(filename,"GE.","")
            name=string.replace(name,"_vs_Others.txt","")
            head=0
            for line in open(InputFile,'rU').xreadlines():
                q=line.rstrip('\r\n')
                q1=string.split(q,'\t')
                if head==0:
                    symbol=q1.index('Symbol')
                    adjp=q1.index('adjp')
                    head=1
                    continue
                else:
                    if q1[symbol] in sfs:
                        mappingdict[name,q1[symbol],"GE"]=q1[adjp]
    dire = export.findParentDir(motifdir)
    output_dir = dire+'MotifResults'
    export.createExportFolder(output_dir)
    output=output_dir+"/Motifresults.txt"
                        
    #output=os.path.join(motifdir,"merged_output_allpvalues_nofold.txt")
    output1=open(output,"w")
    #output1.write("signature"+"\t"+"gene"+"\t"+"tool"+"\t"+"p-value"+"\n")
    for name,gene,key in mappingdict:
        output1.write(name+"\t"+gene+"\t"+key+"\t"+mappingdict[name,gene,key]+"\n")
    output1.close()
    return  output


def Mergeresults(filename):
    Newlist=defaultdict(list)
    Newval={}
    genelst=defaultdict(list)
    Allcomp={}
    dire = export.findParentDir(filename)
    
    output=dire+"/Motifresults_merged.txt"
    #MergedOutput="/Volumes/Pass/MotifAnalyses/Bridger/Exons_MotifAnalyses/merged_output_allpvalues_nofold.txt"
    
    #output="/Volumes/Pass/MotifAnalyses/Bridger/Exons_MotifAnalyses/merged_output_allpvalues_nofold_upd.txt"
    output1=open(output,"w")
    
    output=dire+"/Motifresults_zscores.txt"
    output2=open(output,"w")
    
    output1.write("signature"+"\t"+"gene"+"\t"+"technique"+"\t"+"p-value"+"\t"+"log-transformed"+"\t"+"signature"+"\t"+"gene"+"\t"+"technique"+"\t"+"p-value"+"\t"+"log-transformed"+"\t"+"signature"+"\t"+"gene"+"\t"+"technique"+"\t"+"p-value"+"\t"+"log-transformed"+"\t"+"signature"+"\t"+"gene"+"\t"+"technique"+"\t"+"p-value"+"\t"+"log-transformed"+"\n")
    output2.write("signature"+"\t"+"gene"+"\t"+"cisbp-zscore"+"\t"+"CLIPseq-zscore"+"\t"+"GE-zscore"+"\n")
    for lin in open(filename,'rU').xreadlines():
        genes=[]
        s=lin.rstrip('\r\n')
 
        s1=string.split(s,'\t')
        sig=s1[0]
        
        if s1[2]=="GE":
            genes=[s1[1]]
            
            
        else:
            genes=string.split(s1[1],":")
        
        tool=s1[2]
        
        if 'Cisbp_denovo' in tool:
            tool="Cisbp_denovo"
        if "UpstreamIntron_known" in sig:
            sig = string.replace(sig,"UpstreamIntron_known","Upstream")
          
        if "Intron_known" in s1[0]:
            sig=string.replace(sig,"Intron_known","Combined_intron_new")
            
        if "Exons_known" in s1[0]:
            sig=string.replace(sig,"Exons_known","Exon")
    
        if "DownstreamIntron_known" in s1[0]:
            sig=string.replace(sig,"DownstreamIntron_known","Downstream")
    
            
        for i in range(len(genes)):
            if tool not in genelst[sig,genes[i].upper()]:
                genelst[sig,genes[i].upper()].append(tool)
                
            Newval[sig,tool,genes[i].upper()]=float(s1[3])
            if tool=="GE":
                sig1="Exon:"+sig
                Newval[sig1,tool,genes[i].upper()]=float(s1[3])
                
                genelst[sig1,genes[i].upper()].append(tool)
                sig1="Combined_intron_new:"+sig
                Newval[sig1,tool,genes[i].upper()]=float(s1[3])
                genelst[sig1,genes[i].upper()].append(tool)
        
    zscoredt={}
    cisbp=[]
    clipseq=[]
    ge=[]
  
    for sig,genes in genelst:
       
        tools=[]
        cisbpact=True
        cisbpden=True
        tools=genelst[sig,genes]
       # if genes=="MBNL1":
       # print tools,sig
        a=len(tools)
        if 'Cisbp_Actual' in tools and 'Cisbp_denovo' in tools:
            a=a-1
            if Newval[sig,"Cisbp_Actual",genes]<Newval[sig,"Cisbp_denovo",genes]:
                cisbpden=False
            else:
                cisbpact=False
                
        pval=0.0
        count=0
        if a>1:
            pval=0.0
            count=0
            if "Cisbp_Actual" in tools and cisbpact:
                
                count+=1
               # print str(Newval[sig,"Cisbp_Actual",genes])
                pval=0.0-math.log10(Newval[sig,"Cisbp_Actual",genes])
                output1.write(sig+"\t"+genes+"\t"+"Cisbp_Actual"+"\t"+str(Newval[sig,"Cisbp_Actual",genes])+"\t"+str(pval)+"\t")
                zscoredt[sig,genes]=[pval,]
                cisbp.append(pval)
            else:
                
                output1.write(sig+"\t"+genes+"\t"+"Cisbp_Actual"+"\t"+"NA"+"\t"+"NA"+"\t")
            if 'Cisbp_denovo' in tools and cisbpden:
                count+=1
                #print str(Newval[sig,"Cisbp_denovo",genes])
                pval=0.0-math.log10(Newval[sig,"Cisbp_denovo",genes])
                output1.write(sig+"\t"+genes+"\t"+"Cisbp_denovo"+"\t"+str(Newval[sig,"Cisbp_denovo",genes])+"\t"+str(pval)+"\t")
                zscoredt[sig,genes]=[pval,]
                cisbp.append(pval)
            else:
                output1.write(sig+"\t"+genes+"\t"+"Cisbp_denovo"+"\t"+"NA"+"\t"+"NA"+"\t")
                
            if (sig,genes) not in zscoredt:
                zscoredt[sig,genes]=[0.0,]
                cisbp.append(0.0)
                
            if "Clipseq" in tools:
                count+=1
                #print str(Newval[sig,"Clipseq",genes])
                pval=0.0-math.log10(Newval[sig,"Clipseq",genes])
                output1.write(sig+"\t"+genes+"\t"+"Clipseq"+"\t"+str(Newval[sig,"Clipseq",genes])+"\t"+str(pval)+"\t")
                zscoredt[sig,genes].append(pval)
                clipseq.append(pval)
            else:
               output1.write(sig+"\t"+genes+"\t"+"Clipseq"+"\t"+"NA"+"\t"+"NA"+"\t")
               zscoredt[sig,genes].append(0.0)
               clipseq.append(0.0)
            if "GE" in tools:
                count+=1
                #print str(Newval[sig,"GE",genes])
                pval=0.0-math.log10(Newval[sig,"GE",genes])
                output1.write(sig+"\t"+genes+"\t"+"GE"+"\t"+str(Newval[sig,"GE",genes])+"\t"+str(pval)+"\n")
                zscoredt[sig,genes].append(pval)
                ge.append(pval)
            else:
                output1.write(sig+"\t"+genes+"\t"+"GE"+"\t"+"NA"+"\t"+"NA"+"\n")
                zscoredt[sig,genes].append(0.0)
                ge.append(0.0)
    meancis=np.mean(cisbp)
    meanclip=np.mean(clipseq)
    meange=np.mean(ge)
    sdcis=np.std(cisbp)
    sdclip=np.std(clipseq)
    sdge=np.std(ge)

    
    for sig,genes in zscoredt:
        scores=[]
        scores=zscoredt[sig,genes]
        if len(scores)==3:
            val1=(float(scores[0])-float(meancis))/float(sdcis)
            val2=(float(scores[1])-float(meanclip))/float(sdclip)
            val3=(float(scores[2])-float(meange))/float(sdge)
            output2.write(sig+"\t"+genes+"\t"+str(val1)+"\t"+str(val2)+"\t"+str(val3)+"\n")
        else:
            print "error in zscore calculation"
            print sig,genes
            
        
    
            #pval=pval/float(count)
            #output1.write(str(pval)+"\n")


    

if __name__ == '__main__':

    import getopt
  
    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Insufficient command line flags supplied."
        sys.exit()
    else:
        analysisType = []

        options, remainder = getopt.getopt(sys.argv[1:],'', ['motifAnalysisdir=','GEadjpdir=','SFlist='])
        for opt, arg in options:
            if opt == '--motifAnalysisdir': motifdir=arg
            elif opt == '--GEadjpdir':GEdir=arg
            elif opt == '--SFlist':Sflist=arg
        
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
    
    filename=parseResultfolders(motifdir,GEdir,Sflist)
    
    Mergeresults(filename)