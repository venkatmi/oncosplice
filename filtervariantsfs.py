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
    genelist=[]
    text_file = open("/Users/meenakshi/Documents/InterestingGenes.txt", "rU")
    lines = text_file.readlines()
    for line in lines:
        lin=line.rstrip('\r\n')
        genelist.append(lin)
    text_file.close()
    text_file1 = open("/Users/meenakshi/Documents/SRSF2_like/SRSF2P95H_SRSF2like_U2AF1like.txt","rU")
    lines = text_file1.readlines()
    text_file1.close()
    
    export=open('/Users/meenakshi/Documents/SRSF2_like/SRSF2P95H_SRSF2like_U2AF1like_interestinggenes.txt','w')
    head=0
    for line in lines:
            if head==0:
                export.write(line)
                head=1
                continue
            
            lin=line.rstrip('\r\n')
            lin=string.split(lin,"\t")
            for i in genelist:
                if i in lin[0]:
                    export.write(line)
                    continue
            
            #export.write("\n")
                
initiaterec()
