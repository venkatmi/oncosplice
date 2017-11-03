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


#reading the splicing factor expression file
SFexp_file = "/Users/meenakshi/Documents/ICGS_guidefiles/corrmatrix_round2_imp.txt"


    
    #writing the splicing factor results
export=open('/Users/meenakshi/Documents/ICGS_guidefiles/corrmatrix_round2_events.txt','w')
    
    #processing the splicing events file
head=True
counter=0

for exp in open(SFexp_file,"rU").xreadlines():
        sfvalues=[]
        exp=exp.rstrip(os.linesep)
        if head:
            export.write('Splicing Factor'+'\n')
            
            
            head = False
            continue
        s = string.split(exp,'\t')
        sfname=s[0]
        count=0
        try:
            sfvalues = map(lambda x: float(x), s[1:])
        except Exception:
            
            for value in s[1:]:
                    try: sfvalues.append(float(value))
                    except Exception:
                        #replace NA with 0
                        sfvalues.append(float(0))
        for i in sfvalues:
            if i >0.4:
                break
        else :
            export.write(sfname+'\n')
        
        