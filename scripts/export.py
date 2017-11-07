###export
#Copyright 2005-2008 J. David Gladstone Institutes, San Francisco California
#Author Nathan Salomonis - nsalomonis@gmail.com

#Permission is hereby granted, free of charge, to any person obtaining a copy 
#of this software and associated documentation files (the "Software"), to deal 
#in the Software without restriction, including without limitation the rights 
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
#copies of the Software, and to permit persons to whom the Software is furnished 
#to do so, subject to the following conditions:

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""This module contains generic methods for creating export paths that include new
multiple nested folders and deleting all files within a directory."""

import os
import sys
import string
import unique
import shutil
import UI
if sys.platform == "win32":
    mode = 'wb' ### writes as in a binary mode versus text mode which introduces extra cariage returns in Windows
else:
    mode = 'w' ### writes text mode which can introduce extra carriage return characters (e.g., /r)
    
def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir)
    return dir_list

def getParentDir(filename):
    filename = string.replace(filename,'//','/')
    filename = string.replace(filename,'\\','/')
    return string.split(filename,'/')[-2]
    
def findParentDir(filename):
    ### :: reverses string
    filename = string.replace(filename,'//','/')
    filename = string.replace(filename,'//','/') ### If /// present
    filename = string.replace(filename,'\\','/')
    x = string.find(filename[::-1],'/')*-1 ### get just the parent directory
    return filename[:x]

def findFilename(filename):
    filename = string.replace(filename,'//','/')
    filename = string.replace(filename,'//','/') ### If /// present
    filename = string.replace(filename,'\\','/')
    x = string.find(filename[::-1],'/')*-1 ### get just the parent directory
    return filename[x:]

def ExportFile(filename):
    filename = string.replace(filename,'//','/')
    filename = string.replace(filename,'\\','/')
    dir = findParentDir(filename)
    try: file_var = createExportFile(filename,dir)
    except RuntimeError:
        isFileOpen(filename,dir)
        file_var = createExportFile(filename,dir)
    return file_var

def customFileMove(old_fn,new_fn):
    old_fn = filepath(old_fn)
    new_fn = filepath(new_fn)
    raw = ExportFile(new_fn)
    for line in open(old_fn,'rU').xreadlines():
        if line[0]!='#': ### Applies to Affymetrix APT data (screws up combat)
            raw.write(line)
    raw.close()
    os.remove(old_fn)

def customFileCopy(old_fn,new_fn):
    old_fn = filepath(old_fn)
    new_fn = filepath(new_fn)
    raw = ExportFile(new_fn)
    for line in open(old_fn,'rU').xreadlines():
        if line[0]!='#': ### Applies to Affymetrix APT data (screws up combat)
            raw.write(line)
    raw.close()
    
def isFileOpen(new_file,dir):
    try: 
        file_open = 'yes'
        dir_list = read_directory(dir)
        if len(dir_list)>0:
            while file_open == 'yes':
                file_open = 'no'
                for file in dir_list:
                    if file in new_file: ###Thus the file is open
                        try:
                            fn=filepath(new_file);
                            file_var = open(fn,mode)
                            """
                            except Exception:
                                try: os.chmod(fn,0777) ### It's rare, but this can be a write issue
                                except Exception:
                                    print "This user account does not have write priveledges to change the file:"
                                    print fn,"\nPlease login as an administrator or re-install the software as a non-admin.";sys.exit()
                                file_var = open(fn,'w')"""
                            file_open = 'no'
                        except IOError:
                            print_out = 'Results file: '+fn+ '\nis open...can not re-write.\nPlease close file and select "OK".'
                            try: UI.WarningWindow(print_out,' OK ');
                            except Exception:
                                print print_out; print 'Please correct (hit return to continue)'
                                inp = sys.stdin.readline()
                            file_open = 'yes'
    except OSError: null = []
            
def createExportFile(new_file,dir):
    try:
        #isFileOpen(new_file,dir)  ###Creates problems on Mac - not clear why
        fn=filepath(new_file)
        file_var = open(fn,mode)
        """except Exception:
            try: os.chmod(fn,0777) ### It's rare, but this can be a write issue
            except Exception:
                print "This user account does not have write priveledges to change the file:"
                print fn,"\nPlease login as an administrator or re-install the software as a non-admin.";sys.exit()
            file_var = open(fn,'w')"""
    except Exception:
        createExportDir(new_file,dir) ###Occurs if the parent directory is also missing
        fn=filepath(new_file)
        file_var = open(fn,mode)
        """except Exception:
            try: os.chmod(fn,0777) ### It's rare, but this can be a write issue
            except Exception:
                print "This user account does not have write priveledges to change the file:"
                print fn,"\nPlease login as an administrator or re-install the software as a non-admin.";sys.exit()
            file_var = open(fn,'w')"""
    return file_var

def createExportDirAlt(new_file,dir):
    ### Original method for creating a directory path that is not present
    ### Works by going backwards (not ideal)
    dir = string.replace(dir,'//','/')
    dir = string.replace(dir,'\\','/')
    dir = string.replace(dir,'\\','/')
    dir_ls = string.split(dir,'/')
    if len(dir_ls) != 1:
        index = 1
        while index < (len(dir_ls)+1):
            parent_dir = string.join(dir_ls[:index],'/')
            index+=1
            try:pfn = filepath(parent_dir); os.mkdir(pfn)  
            except OSError: continue
        createExportFile(new_file,dir)
    else: print "Parent directory not found locally for", dir_ls

def createExportDir(new_file,dir):
    ### New method for creating a directory path that is not present
    ### Works by going forward (check if a base path is present and then go up)
    dir = string.replace(dir,'//','/')
    dir = string.replace(dir,'\\','/')
    dir = string.replace(dir,'\\','/')
    dir = filepath(dir)
    dir_ls = string.split(dir,'/')
    i = 1; paths_added = 'no'
    while i <= len(dir_ls):
        new_dir = string.join(dir_ls[:i],'/')
        status = verifyDirectory(new_dir)
        if status == 'no':
            try: os.mkdir(new_dir); paths_added = 'yes'
            except Exception: paths_added = 'yes'
        i+=1
    if paths_added == 'yes':
        try:
            fn=filepath(new_file)
            file_var = open(fn,mode)
        except Exception:
            print "Parent directory not found locally for", [dir,new_file]
    #else: print "Parent directory not found locally for", [dir,new_file]; sys.exit()
    
def createDirPath(dir):
    ### New method for creating a directory path that is not present
    ### Works by going forward (check if a base path is present and then go up)
    dir = string.replace(dir,'//','/')
    dir = string.replace(dir,'\\','/')
    dir = string.replace(dir,'\\','/')
    dir = filepath(dir)
    dir_ls = string.split(dir,'/')
    i = 1; paths_added = 'no'
    while i <= len(dir_ls):
        new_dir = string.join(dir_ls[:i],'/')
        status = verifyDirectory(new_dir)
        if status == 'no':
            try: os.mkdir(new_dir); paths_added = 'yes'
            except Exception: paths_added = 'yes'
        i+=1

def verifyDirectory(dir):
    try: dir_list = read_directory(dir); verified = 'yes'
    except Exception: verified = 'no'
    #print 'verify',[dir],verified
    if verified == 'no': ### Can occur for the AltDatabase dir, since EnsMart62 will be added an not found
        try: dir_list = os.listdir(dir); verified = 'yes'
        except Exception: verified = 'no'
    return verified
    
def createExportFolder(dir):
    dir = string.replace(dir,'//','/')
    dir = string.replace(dir,'\\','/')
    dir = string.replace(dir,'\\','/')
    dir_ls = string.split(dir,'/')
    if len(dir_ls) != 1:
        index = 1
        while index < (len(dir_ls)+1):
            parent_dir = string.join(dir_ls[:index],'/')
            index+=1
            #print "Trying to create the directory:",parent_dir
            try:
                pfn = filepath(parent_dir)
                #print pfn
                os.mkdir(pfn)
                #print parent_dir, 'written' #, index-1, len(dir_ls)
            except OSError:
                #print "Can not write this dir"
                #break
                continue
    else: print "Parent directory not found locally for", dir_ls

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def cleanFile(source_file,removeExtra=None):
    ### Some files have extra odd encoding that results in blank new lines in the extracted file
    ### For succeptible directories copy all files line by line, removing existing end of lines
    file = findFilename(source_file); temp_file = 'tempdir/'+file
    data = ExportFile(temp_file)
    fn=filepath(source_file)
    for line in open(fn,'rU').xreadlines():
        line = cleanUpLine(line)
        writeFile=True
        if removeExtra!=None:
            if line[0]==removeExtra: writeFile=False
        if len(line)>0 and writeFile: data.write(line+'\n')
    data.close()
    ### Replace old file with new file
    copyFile(temp_file,source_file)
    os.remove(temp_file)
    print 'copied',file
    
def copyFile(source_file,destination_file):
    dir = findParentDir(destination_file)
    try: createExportFolder(dir)
    except Exception: null=[] ### already exists
    shutil.copyfile(source_file,destination_file)
    #print '\nFile copied to:',destination_file
    
def deleteFolder(dir):
    try:
        dir = filepath(dir); dir_list = read_directory(dir) ### Get all files in directory
        #try: print 'deleting dir:',dir
        #except Exception: null=None ### Occurs due to a Tkinter issue sometimes
        for file in dir_list:
            fn = filepath(dir+'/'+file)
            try:
                if '.' in file: os.remove(fn)
                else: deleteFolder(fn) ### Remove subdirectories
            except Exception: None
        os.removedirs(dir)
        #print dir
        return 'success'
    except OSError: return 'failed'

if __name__ == '__main__':
    customFileCopy('/Volumes/SEQ-DATA/IlluminaBodyMap/pooled/pooled/Thyroid-ERR030872__exons.bed','/Volumes/SEQ-DATA/IlluminaBodyMap/pooled/Thyroid-ERR030872__exons.bed');sys.exit()
    createExportFolder('Databases/null/null'); sys.exit()
    createExportDir('C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/RNASeq/hESC-NP/TopHat-hESC_differentiation/AltExpression/pre-filtered/counts/a.txt','C:/Users/Nathan Salomonis/Desktop/Gladstone/1-datasets/RNASeq/hESC-NP/TopHat-hESC_differentiation/AltExpression/pre-filtered/counts'); sys.exit()
    deleteFolder('BuildDBs/Entrez/Gene2GO');kill
    fn = '/Users/nsalomonis/Desktop/GSE13297_RAW//AltExpression/ExonArray/Hs/Hs_Exon_CS_vs_hESC.p5_average.txt'
    createExportFile(fn,'/Users/nsalomonis/Desktop/GSE13297_RAW//AltExpression/ExonArray/Hs')
    ExportFile(fn)
