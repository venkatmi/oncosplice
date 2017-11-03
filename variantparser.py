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
    text_file = open("/Volumes/MyPassport/Leucegene_cosmic_tcga_exactamtch/Cosmic_moderate_result-upd.txt", "rU")
    lines = text_file.readlines()
    
    text_file.close()
    export=open('/Volumes/MyPassport/Leucegene_cosmic_tcga_exactamtch/Cosmic_moderate_result-upd_parsed.txt','w')
    head=0
    for line in lines:
            if head==0:
                export.write(line)
                head=1
                continue
            newlin=[]
            lin=line.rstrip('\r\n')
            lin=string.split(lin,"\t")
            export.write(lin[0]+"\t")
            
            newlin=lin[1:15]
            #print newlin;sys.exit()
            for i in lin[15:]:
                
                if i=="./.":
                    newlin.append(str(0))
                else:
                    newlin.append(str(1))
            
            export.write("\t".join(newlin))
            export.write("\n")
                
initiaterec()
