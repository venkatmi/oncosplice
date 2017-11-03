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
head=0
import math


def initiaterec():
    
    text_file = open("/Volumes/MyPassport/Leucegene/Leucegene_Results/AlternativeOutput/ExpressionInput/exp.LeucegeneSplicing.txt", "rU")
    lines = text_file.readlines()
    
    text_file.close()
    export=open('/Volumes/MyPassport/Leucegene/Leucegene_Results/AlternativeOutput/ExpressionInput/exp.LeucegeneSplicing_v32std.txt','w')
    head=0
    for line in lines:
        list1=[];list2=[];
        if head ==0:
            head=1
            export.write(line)
            continue
        else:
            lin=line.rstrip('\r\n')
            t = string.split(lin,'\t')
            export.write(t[0])
            for i in range(1,len(t)):
                if t[i]=="" or t[i]=="0":
                    continue
                else:
                    list1.append(float(t[i]))
            std_val=numpy.std(list1)
            std_val2=float(2.0*std_val)
            for i in range(1,len(t)):
                if t[i]=="":
                    list2.append('0.0')
                else:
                    if float(t[i]) > std_val2 or float(t[i]) <-std_val2:
                        list2.append(round(float(t[i]),2))
                    else:
                        list2.append('0.0')
                
            for i in list2:
                export.write("\t"+str(i))
            export.write('\n')
                    
                
initiaterec()

