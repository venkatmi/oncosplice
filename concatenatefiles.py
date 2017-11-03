#!/usr/bin/env python

import os

def get_filepaths(directory):
    """
    This function will generate the file names in a directory 
    tree by walking the tree either top-down or bottom-up. For each 
    directory in the tree rooted at directory top (including top itself), 
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        output=("output_exp.txt")
        export=open(output,"w")
        for filename in files:
            # Join the two strings in order to form the full filepath.
         if "DS" not in filename:
            filepath = os.path.join(root, filename)
            lst=[]
           # file_paths.append(filepath)  # Add it to the list.
            for line in open(filepath,'rU').xreadlines():
                line = line.rstrip('\n')
                lst.append(line)
            
            
            export.write(filename+'\t')
            export.write("\t".join(lst))
            export.write('\n')
 # Self-explanatory.

get_filepaths("/Users/meenakshi/Documents/KNN_classification/GeneExpression/results/SVC_proba")