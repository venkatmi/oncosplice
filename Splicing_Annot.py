#!/usr/bin/env python
import numpy as np
import pylab as pl
import sys,string
import os
from numpy import loadtxt
from collections import Counter
from collections import defaultdict

#Annotate PSI file's splicing events using alternate_juntion  and alternate_junction_de-novo

def parse_junctionfiles():
    head=0;head1=0;
    count=0;count1=0;
    psievents={}
    # Create a dictionary for the psi event dict[min_isoform,major_isoform] and value will be the annotation that will assigned later in the code
    for line in open('/Volumes/MyPassport/Leucegene/Leucegene_Results/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI.txt','rU').xreadlines():
        line = line.rstrip('\n')
        line=string.replace(line,"_",".")
        flag=0
        if head ==0:
            head=1
            continue
        psi=string.split(line,'\t')
        
        if (psi[2],psi[3] not in psievents) and (psi[3],psi[2] not in psievents):
            psievents[psi[2],psi[3]]=''
            psievents[psi[3],psi[2]]=''
    print"psi events dictionary created"
    
    # Parse to the de-novo junction file and update the value of the psi events in the dictionary 
    head=0   
    for line in open('/Volumes/MyPassport/Leucegene_bedfiles/AltDatabase/EnsMart72/ensembl/Hs/Hs_alternative_junctions_de-novo.txt','rU').xreadlines():
        line = line.rstrip('\n')
        flag=0
        if head ==0:
           head=1
           continue
        de_event=string.split(line,'\t')
        min_isoform=de_event[0]+':'+de_event[2]
        maj_isoform=de_event[0]+":"+de_event[3]
       
        if (min_isoform,maj_isoform in psievents) or (maj_isoform,min_isoform in psievents):
            
            psievents[min_isoform,maj_isoform]=de_event[4]
            psievents[maj_isoform,min_isoform]=de_event[4]
    print "denovo junction file parsed"
    
    # Parse to the junction file and update the value of the psi events in the dictionary 
    head=0  
    for line in open('/Users/meenakshi/Downloads/AltAnalyze_v.2.0.9.4-Py/AltDatabase/EnsMart72/ensembl/Hs/Hs_alternative_junctions.txt','rU'):
        line = line.rstrip('\n')
        flag=0
        if head ==0:
           head=1
           continue
        ju_event=string.split(line,'\t')
        min_isoform=ju_event[0]+':'+ju_event[2]
        maj_isoform=ju_event[0]+":"+ju_event[3]
        if (min_isoform,maj_isoform in psievents) or (maj_isoform,min_isoform in psievents):
            psievents[min_isoform,maj_isoform]=ju_event[4]
            psievents[maj_isoform,min_isoform]=ju_event[4]
    print "main junction file parsed"
    
    # write the updated psi file with the annotations into a new file and annotate events that have not been annotated by the junction files
    head=0
    export=open('/Volumes/MyPassport/Leucegene/Leucegene_Results/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSIupd1.txt','w')
    count=0
    for line in open('/Volumes/MyPassport/Leucegene/Leucegene_Results/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI.txt','rU').xreadlines():
        line = line.rstrip('\n')
        flag=0
        if head ==0:
            export.write(line+'\t'+'splicingevents'+'\n')
            head=1
            continue
        psiEvent=string.split(line,'\t')
        psiEvent_min=string.replace(psiEvent[2],"_",".")
        psiEvent_maj=string.replace(psiEvent[3],"_",".")
        try:
            if psievents[psiEvent_min,psiEvent_maj]=='':
                min_exons=string.split(psiEvent[2],":")
                maj_exons=string.split(psiEvent[3],":")
                if len(min_exons)>2 or len(maj_exons)>2:
                    export.write(line+'\t'+'trans-splicing'+'\n')
                    continue
                else:
                    min_exonspos=string.split(min_exons[1],"-")
                    maj_exonspos=string.split(maj_exons[1],"-")
                    if ('I' in min_exons[1]) or ('I' in maj_exons[1]):
                        if ('I1.' in min_exonspos[0]) or ('I1.' in maj_exonspos[0]):
                            export.write(line+'\t'+'alt-N-term'+'\n')
                            continue
                        else:
                            export.write(line+'\t'+'cassette-exon'+'\n')
                            continue
                    min_exonpos1=string.split(min_exonspos[0],".")[0]
                    min_exonpos2=string.split(min_exonspos[1],".")[0]
                    maj_exonpos1=string.split(maj_exonspos[0],".")[0]
                    maj_exonpos2=string.split(maj_exonspos[1],".")[0]
                    if min_exonpos1==min_exonpos2 or maj_exonpos1==maj_exonpos2:
                            export.write(line+'\t'+'exon-region-exclusion'+'\n')
                            continue
                export.write(line+'\t'+''+'\n')
            else:
                export.write(line+'\t'+psievents[psiEvent_min,psiEvent_maj]+'\n')
        except Exception():
            export.write(line+'\t'+''+'\n')


parse_junctionfiles()
