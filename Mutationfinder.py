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
mutdict = {}
import sys
head=0
import math


def initiaterec(lower,upper,aggregate="false"):
    samplelst=[]
    mutfile=""
    export=open('/Volumes/MyPassport/mutations-list.txt','w')
    head=0
    for exp1 in open(mutfile,"rU").xreadlines():
        lin=exp1.rstrip('\r\n')
        lin=string.split(lin,"\t")
        if head==0:
            export.write(exp1)
            samplelst=lin
            head=1
            continue
        newlin=[]
        
        k=len(lin)-1
        if k >lower and k< upper:
            if i ==0:
                   newlin.append(str(0))
            else:
                    i=-float(i)
                    newlin.append(str(i))
                    
            
            export.write("\t".join(newlin))
            export.write("\n")
initiaterec()

