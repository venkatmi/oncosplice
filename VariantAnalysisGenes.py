#!/usr/bin/env python
import os
import os.path
global test_list
import matplotlib.pyplot as plt
import copy
from collections import Counter
from collections import defaultdict
import string
import numpy
tandem = {}
import sys
head=0
import math


def initiaterec():
    text_file = open("/Volumes/MyPassport/Leucegene/combined_all_samples/Moderate_inputfile-output.txt", "rU")
    lines = text_file.readlines()
    
    text_file.close()
    export=open('/Volumes/MyPassport/Leucegene/combined_all_samples/Moderate_inputfile-output_genes.txt','w')
    head=0
    pat=[]
    patdict=defaultdict(list)
    for line in lines:
        if head==0:
            export.write(line)
            lin=line.rstrip('\r\n')
            lin=string.split(lin,"\t")
            pat=lin[1:]
            head=1
            continue
        newlin=[]
        pastlin=[]
        values=[]
        lin=line.rstrip('\r\n')
        lin=string.split(lin,"\t")
        #export.write(lin[0]+"\t")
        newlin=string.split(lin[0],":")[0]
        newlin=string.split(newlin,"|")
        newlin=set(newlin)
        print newlin
        values=lin[1:]
        for i in newlin:
            if i not in patdict:
                for j in range(len(values)):
                    if values[j]=='1':
                        if j not in patdict[i]:
                            patdict[i].append(j)
            else:
                for j in range(len(values)):
                    if values[j]=='1':
                        if j not in patdict[i]:
                            patdict[i].append(j)
    val=[]          
    for key in patdict:
        val=patdict[key]
        export.write(key)
        for j in range(len(pat)):
            if j in val:
                export.write('\t'+'1')
            else:
                export.write('\t'+'0')
        export.write("\n")
                
initiaterec()