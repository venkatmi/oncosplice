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
    text_file1 = open("/Users/meenakshi/Documents/exp.splicing.txt","rU")
    lines = text_file1.readlines()
    text_file1.close()
    
    export=open('/Volumes/MyPassport/Leuce-upd/exp.splicing-output.txt','w')
    head=0
    for line in lines:
            if head==0:
                export.write(line)
                print line
                head=1
                continue
            newlin=[]
            lin=line.rstrip('\r\n')
            lin=string.split(lin,"\t")
            export.write(lin[0]+"\t")
            for i in lin[1:]:
                if i ==0:
                   newlin.append(str(0))
                else:
                    i=-float(i)
                    newlin.append(str(i))
                
            
            export.write("\t".join(newlin))
            export.write("\n")
initiaterec()