### Based on code from AltAnalyze's LineageProfiler (http://altanalyze.org)
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


"""
This script iterates the LineageProfiler algorithm (correlation based classification method) to identify sample types relative to one
of two references given one or more gene models. The main function is runLineageProfiler.

The program performs the following actions:
    1) Import a tab-delimited reference expression file with three columns (ID, biological group 1, group 2) and a header row (biological group names)
    2) Import a tab-delimited expression file with gene IDs (column 1), sample names (row 1) and normalized expression values (e.g., delta CT values)
    3) (optional - import existing models) Import a tab-delimited file with comma delimited gene-models for analysis
    4) (optional - find new models) Identify all possible combinations of gene models for a supplied model size variable (e.g., --s 7)
    5) Iterate through any supplied or identified gene models to obtain predictions for novel or known sample types
    6) Export prediction results for all analyzed models to the folder SampleClassification.
    7) (optional) Print the top 20 scores and models for all possible model combinations of size --s
    
"""
    
import sys, string
import math
import os.path
import copy
import time
import getopt
import scipy

try: import unique ### Not required (used in AltAnalyze)
except Exception: None
try: import export ### Not required (used in AltAnalyze)
except Exception: None

#import salstat_stats; reload(salstat_stats)
try:
    from scipy import stats
    use_scipy = True
except Exception:
    use_scipy = False ### scipy is not required but is used as a faster implementation of Fisher Exact Test when present

def filepath(filename):
    try: fn = unique.filepath(filename)
    except Exception: fn = filename
    return fn

def exportFile(filename):
    try: export_data = export.ExportFile(filename)
    except Exception: export_data = open(filename,'w')
    return export_data

def makeUnique(item):
    db1={}; list1=[]; k=0
    for i in item:
        try: db1[i]=[]
        except TypeError: db1[tuple(i)]=[]; k=1
    for i in db1:
        if k==0: list1.append(i)
        else: list1.append(list(i))
    list1.sort()
    return list1

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def returnLargeGlobalVars():
    ### Prints all large global variables retained in memory (taking up space)
    all = [var for var in globals() if (var[:2], var[-2:]) != ("__", "__")]
    for var in all:
        try:
            if len(globals()[var])>1:
                print var, len(globals()[var])
        except Exception: null=[]
        
def clearObjectsFromMemory(db_to_clear):
    db_keys={}
    for key in db_to_clear: db_keys[key]=[]
    for key in db_keys:
        try: del db_to_clear[key]
        except Exception: 
            try:
                for i in key: del i ### For lists of tuples
            except Exception: del key ### For plain lists

def int_check(value):
    val_float = float(value)
    val_int = int(value)
    if val_float == val_int:
        integer_check = 'yes'
    if val_float != val_int:
        integer_check = 'no'
    return integer_check

def IQR(array):
    k1 = 75
    k2 = 25
    array.sort()
    n = len(array)
    value1 = float((n*k1)/100)
    value2 = float((n*k2)/100)
    if int_check(value1) == 'no':
        k1_val = int(value1) + 1
    if int_check(value1) == 'yes':
        k1_val = int(value1)
    if int_check(value2) == 'no':
        k2_val = int(value2) + 1
    if int_check(value2) == 'yes':
        k2_val = int(value2)
    median_val = scipy.median(array)
    upper75th = array[k1_val]
    lower25th = array[k2_val]
    int_qrt_range = upper75th - lower25th
    T1 = lower25th-(1.5*int_qrt_range)
    T2 = upper75th+(1.5*int_qrt_range)
    return lower25th,median_val,upper75th,int_qrt_range,T1,T2

class IQRData:
    def __init__(self,maxz,minz,medz,iq1,iq3):
        self.maxz = maxz; self.minz = minz
        self.medz = medz; self.iq1 = iq1
        self.iq3 = iq3
    def Max(self): return self.maxz
    def Min(self): return self.minz
    def Medium(self): return self.medz
    def IQ1(self): return self.iq1
    def IQ3(self): return self.iq3
    def SummaryValues(self):
        vals = string.join([str(self.IQ1()),str(self.Min()),str(self.Medium()),str(self.Max()),str(self.IQ3())],'\t')
        return vals
    
def importGeneModels(geneModels):
    fn=filepath(geneModels); x=0
    geneModels=[]
    for line in open(fn,'rU').xreadlines():
        genes = cleanUpLine(line)
        genes = string.replace(genes,"'",'')
        genes = string.replace(genes,' ',',')
        genes = string.split(genes,',')
        models=[]
        for gene in genes:
            if len(gene)>0:
                models.append(gene)
        if len(models)>0:
            geneModels.append(models)
    return geneModels
                
######### Below code deals is specific to this module #########
def runLineageProfiler(species,array_type,exp_input,exp_output,codingtype,compendium_platform,modelSize=None,customMarkers=False,geneModels=False):
    """ This code differs from LineageProfiler.py in that it is able to iterate through the LineageProfiler functions with distinct geneModels
    that are either supplied by the user or discovered from all possible combinations. """
    
    global exp_output_file; exp_output_file = exp_output; global targetPlatform
    global tissue_specific_db; global expession_subset; global tissues; global sample_headers
    global analysis_type; global coding_type; coding_type = codingtype
    global tissue_to_gene; tissue_to_gene = {}; global platform; global cutoff
    global customMarkerFile; global delim
    customMarkerFile = customMarkers
    if geneModels == False: geneModels = []
    else:
        geneModels = importGeneModels(geneModels)
    
    if '\\' in exp_input: delim = '\\'
    elif '//' in exp_input: delim = '//'
    else: delim = "/"
    
    print '\nRunning LineageProfiler analysis on',string.split(exp_input,delim)[-1][:-4]
    
    global correlate_by_order; correlate_by_order = 'no'
    global rho_threshold; rho_threshold = -1
    global correlate_to_tissue_specific; correlate_to_tissue_specific = 'no'
    platform = array_type
    cutoff = 0.01
    global value_type
    
    if 'stats.' in exp_input:
        value_type = 'calls'
    else:
        value_type = 'expression'
    
    tissue_specific_db={}; expession_subset=[]; sample_headers=[]; tissues=[]
    if len(array_type)==2:
        ### When a user-supplied expression is provided (no ExpressionOutput files provided - importGeneIDTranslations)
        vendor, array_type = array_type
        platform = array_type
    else: vendor = 'Not needed'

    if 'RawSplice' in exp_input or 'FullDatasets' in exp_input or coding_type == 'AltExon':
        analysis_type = 'AltExon'
        if platform != compendium_platform: ### If the input IDs are not Affymetrix Exon 1.0 ST probesets, then translate to the appropriate system
            translate_to_genearray = 'no'
            targetPlatform = compendium_platform
            translation_db = importExonIDTranslations(array_type,species,translate_to_genearray)
            keyed_by = 'translation'
        else: translation_db=[]; keyed_by = 'primaryID'; targetPlatform = compendium_platform
    elif array_type == "3'array" or array_type == 'AltMouse':
        ### Get arrayID to Ensembl associations
        if vendor != 'Not needed':
            ### When no ExpressionOutput files provided (user supplied matrix)
            translation_db = importVendorToEnsemblTranslations(species,vendor,exp_input)
        else:
            translation_db = importGeneIDTranslations(exp_output)
        keyed_by = 'translation'
        targetPlatform = compendium_platform
        analysis_type = 'geneLevel'
    else:
        translation_db=[]; keyed_by = 'primaryID'; targetPlatform = compendium_platform; analysis_type = 'geneLevel'

    targetPlatform = compendium_platform ### Overides above
    try: importTissueSpecificProfiles(species)
    except Exception:
        try:
            try:
                targetPlatform = 'exon'
                importTissueSpecificProfiles(species)
            except Exception:
                try:
                    targetPlatform = 'gene'
                    importTissueSpecificProfiles(species)
                except Exception: 
                    targetPlatform = "3'array"
                    importTissueSpecificProfiles(species)
        except ZeroDivisionError:
            print 'No compatible compendiums present...'
            print e
            forceError
                
    all_marker_genes=[]            
    for gene in tissue_specific_db:
        all_marker_genes.append(gene)

    if len(geneModels)>0:
        allPossibleClassifiers = geneModels
    elif modelSize != None:
        ### A specific model size has been specified (e.g., find all 10-gene models)
        allPossibleClassifiers = getRandomSets(all_marker_genes,modelSize)
    else:
        allPossibleClassifiers = [all_marker_genes]
        

    num=1
    if len(allPossibleClassifiers)<16:
        print 'Using:'
        for model in allPossibleClassifiers:
            print 'model',num,model
            num+=1

    hit_list=[]
    
    print 'Number of references to compare to:',len(tissues), tissues
    ### Iterate through LineageProfiler for all gene models (allPossibleClassifiers)
    times = 1; k=1000; l=1000; hits=[]; fails=[]; f=0; s=0; sample_diff_z={}; prognostic_class1_db={}; prognostic_class2_db={}
    begin_time = time.time()
    for classifiers in allPossibleClassifiers:
        tissue_to_gene={}; expession_subset=[]; sample_headers=[]; classifier_specific_db={}
        for gene in classifiers:
            try: classifier_specific_db[gene] = tissue_specific_db[gene]
            except Exception: None
        importGeneExpressionValues(exp_input,classifier_specific_db,keyed_by,translation_db)
        ### If the incorrect gene system was indicated re-run with generic parameters
        if len(expession_subset)==0 and (array_type == "3'array" or array_type == 'AltMouse'):
            translation_db=[]; keyed_by = 'primaryID'; targetPlatform = compendium_platform; analysis_type = 'geneLevel'
            importGeneExpressionValues(exp_input,tissue_specific_db,keyed_by,translation_db)
        if len(expession_subset)!=len(classifiers): f+=1
        if len(expession_subset)==len(classifiers): ### Sometimes a gene or two are missing from one set
            s+=1
            zscore_output_dir,tissue_scores = analyzeTissueSpecificExpressionPatterns()
            #except Exception: print len(classifier_specific_db), classifiers; error
            headers = list(tissue_scores['headers']); del tissue_scores['headers']
            if times == k:
                end_time = time.time()
                print int(end_time-begin_time),'seconds'
                k+=l
            times+=1; index=0; positive=0; positive_score_diff=0
            sample_number = (len(headers)-1)
            population1_denom=0; population1_pos=0; population2_pos=0; population2_denom=0
            diff_positive=[]; diff_negative=[]
            while index < sample_number:
                diff_z = tissue_scores[tissues[0]][index]-tissue_scores[tissues[-1]][index]
                try: sample_diff_z[headers[index+1]].append(diff_z)
                except Exception: sample_diff_z[headers[index+1]] =[diff_z]
                if headers[index+1] not in prognostic_class1_db:
                    prognostic_class1_db[headers[index+1]]=0 ### Create a default value for each sample
                if headers[index+1] not in prognostic_class2_db:
                    prognostic_class2_db[headers[index+1]]=0 ### Create a default value for each sample
                if diff_z>0:
                    prognostic_class1_db[headers[index+1]]+=1
                if diff_z<0:
                    prognostic_class2_db[headers[index+1]]+=1
                if diff_z>0 and (tissues[0]+':' in headers[index+1]):
                    positive+=1; positive_score_diff+=abs(diff_z)
                    population1_pos+=1; diff_positive.append(abs(diff_z))
                    hits.append(headers[index+1]) ### see which are correctly classified
                elif diff_z<0 and (tissues[-1]+':' in headers[index+1]):
                    positive+=1; positive_score_diff+=abs(diff_z)
                    population2_pos+=1; diff_positive.append(abs(diff_z))
                    hits.append(headers[index+1]) ### see which are correctly classified
                elif diff_z>0 and (tissues[-1]+':' in headers[index+1]): ### Incorrectly classified
                    diff_negative.append(abs(diff_z))
                    fails.append(headers[index+1])
                elif diff_z<0 and (tissues[0]+':' in headers[index+1]): ### Incorrectly classified
                    #print headers[index+1]
                    diff_negative.append(abs(diff_z))
                    fails.append(headers[index+1])
                if (tissues[0]+':' in headers[index+1]):
                    population1_denom+=1
                else:
                    population2_denom+=1
                index+=1
            percent_positive = (float(positive)/float(index))*100
            hit_list.append([percent_positive,population1_pos, population1_denom,population2_pos,population2_denom,[avg(diff_positive),avg(diff_negative)],positive_score_diff,classifiers])

    root_dir = string.join(string.split(exp_output_file,'/')[:-1],'/')+'/'
    dataset_name = string.replace(string.split(exp_input,'/')[-1][:-4],'exp.','')
    output_classification_file = root_dir+'SampleClassification/'+dataset_name+'-SampleClassification.txt'
    try: os.mkdir(root_dir+'SampleClassification')
    except Exception: None
    export_summary = exportFile(output_classification_file)
    models = []
    for i in allPossibleClassifiers:
        i = string.replace(str(i),"'",'')[1:-1]
        models.append(i)
    headers = string.join(['Samples']+[tissues[0]+' Predicted Hits',tissues[1]+' Predicted Hits']+['Overall Prognostic Score','Median Z-score']+models,'\t')+'\n'
    export_summary.write(headers)
    for sample in prognostic_class1_db:
        class1_score = prognostic_class1_db[sample]
        class2_score = prognostic_class2_db[sample]
        zscore_distribution = map(str,sample_diff_z[sample])
        dist_list = map(float,zscore_distribution) ### convert to list
        median_score = scipy.median(dist_list)
        values = [sample,str(class1_score),str(class2_score),str(class1_score-class2_score),str(median_score)]
        #print values
        #print zscore_distribution
        export_summary.write(string.join(values+zscore_distribution,'\t')+'\n')
    export_summary.close()
    print 'Results file written to:',root_dir+'SampleClassification/'+dataset_name+'-SampleClassification.txt','\n'

    hit_list.sort(); hit_list.reverse()
    top_hit_list=[]
    top_hit_db={}
    hits_db={}; fails_db={}
    
    for i in sample_diff_z:
        zscore_distribution = sample_diff_z[i]
        maxz = max(zscore_distribution); minz = min(zscore_distribution)
        lower25th,medz,upper75th,int_qrt_range,T1,T2 = IQR(zscore_distribution)
        if float(maxz)>float(T2): maxz = T2
        if float(minz) < float(T1): minz = T1
        iqr = IQRData(maxz,minz,medz,lower25th,upper75th)
        #sample_diff_z[i] = iqr
        sample_diff_z[i] = string.join(map(str,zscore_distribution),'\t')
        
    for i in hits:
        try: hits_db[i]+=1
        except Exception: hits_db[i]=1
    for i in fails:
        try: fails_db[i]+=1
        except Exception: fails_db[i]=1
    for i in fails_db:
        if i not in hits:
            try:
                #print i+'\t'+'0\t'+str(fails_db[i])+'\t'+ sample_diff_z[i]
                None
            except Exception:
                #print i
                None

    if modelSize != None:
        hits=[]
        for i in hits_db:
            hits.append([hits_db[i],i])
        hits.sort()
        hits.reverse()
        for i in hits:
            if i[1] in fails_db: fail = fails_db[i[1]]
            else: fail = 0
            try:
                #print i[1]+'\t'+str(i[0])+'\t'+str(fail)+'\t'+sample_diff_z[i[1]]
                None
            except Exception:
                #print i[1]
                None
        
        for i in hit_list:
            if i[0]>0:
                top_hit_list.append(i[-1])
                top_hit_db[i[-1]]=i[0]
                
    print 'Top 20 top hits'
    for i in hit_list[:20]:
        print i[:5],i[-1]

    return top_hit_db

def factorial(n):
    ### Code from http://docs.python.org/lib/module-doctest.html
    if not n >= 0:
        raise ValueError("n must be >= 0")
    if math.floor(n) != n:
        raise ValueError("n must be exact integer")
    if n+1 == n:  # catch a value like 1e300
        raise OverflowError("n too large")
    result = 1
    factor = 2
    while factor <= n:
        result *= factor
        factor += 1
    return result

def choose(n,x):
    """Equation represents the number of ways in which x objects can be selected from a total of n objects without regard to order."""
    #(n x) = n!/(x!(n-x)!)
    f = factorial
    result = f(n)/(f(x)*f(n-x))
    return result

def getRandomSets(a,size):
    #a = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    #size = 4
    select_set={'ENSG00000140678':'ITGAX','ENSG00000105835':'NAMPT','ENSG00000027697':'IFNGR1','ENSG00000120129':'DUSP1','ENSG00000003402':'CFLAR','ENSG00000113269':'RNF130'}
    select_set={}
    
    select_set2={'ENSG00000163602': 'RYBP'}
    negative_select = {'ENSG00000105352':'CEACAM4'}
    negative_select={}

    import random
    possible_sets = choose(len(a),size)
    print 'Possible',size,'gene combinations to test',possible_sets,
    permute_ls = []; done = 0; permute_db={}
    while done == 0:
        b = list(tuple(a)); random.shuffle(b)
        bx_set={}
        i = 0
        while i < len(b):
            try:
                bx = b[i:i+size]; bx.sort()
                if len(bx)==size: permute_db[tuple(bx)]=None
                else: break
            except Exception: break
            i+=1
        if len(permute_db) == possible_sets:
            done=1; break
    for i in permute_db:
        add=0; required=0; exclude=0
        for l in i:
            if len(select_set)>0:
                if l in select_set: add+=1
                #if l in select_set2: required+=1
            #if l in negative_select: exclude+=1
            else: add = 1000
        if add>2 and exclude==0:# and required==1:
            permute_ls.append(i)
    #print len(permute_ls)
    return permute_ls

def importVendorToEnsemblTranslations(species,vendor,exp_input):
    translation_db={}
    """
    ### Faster method but possibly not as good
    uid_db = simpleUIDImport(exp_input)
    import gene_associations
    ### Use the same annotation method that is used to create the ExpressionOutput annotations
    array_to_ens = gene_associations.filterGeneToUID(species,'Ensembl',vendor,associated_IDs)
    for arrayid in array_to_ens:
        ensembl_list = array_to_ens[arrayid]
        try: translation_db[arrayid] = ensembl_list[0] ### This first Ensembl is ranked as the most likely valid based on various metrics in getArrayAnnotationsFromGOElite
        except Exception: None
    """
    translation_db={}
    import BuildAffymetrixAssociations
    
    ### Use the same annotation method that is used to create the ExpressionOutput annotations
    use_go = 'yes'
    conventional_array_db={}
    conventional_array_db = BuildAffymetrixAssociations.getArrayAnnotationsFromGOElite(conventional_array_db,species,vendor,use_go)
    for arrayid in conventional_array_db:
        ca = conventional_array_db[arrayid]
        ens = ca.Ensembl()
        try: translation_db[arrayid] = ens[0] ### This first Ensembl is ranked as the most likely valid based on various metrics in getArrayAnnotationsFromGOElite
        except Exception: None
    
    return translation_db

def importTissueSpecificProfiles(species):
    if analysis_type == 'AltExon':
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_'+targetPlatform +'_tissue-specific_AltExon_protein_coding.txt'
    else:
        filename = 'AltDatabase/ensembl/'+species+'/'+species+'_'+targetPlatform +'_tissue-specific_'+coding_type+'.txt'
    if customMarkerFile != False:
        filename = customMarkerFile

    if value_type == 'calls':
        filename = string.replace(filename,'.txt','_stats.txt')
    fn=filepath(filename); x=0
    tissues_added={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        
        if x==0:
            print 'Importing the tissue compedium database:',string.split(filename,delim)[-1][:-4]
            headers = t; x=1; index=0
            for i in headers:
                if 'UID' == i: ens_index = index
                if analysis_type == 'AltExon': ens_index = ens_index ### Assigned above when analyzing probesets
                elif 'Ensembl' in i: ens_index = index
                if 'marker-in' in i: tissue_index = index+1; marker_in = index
                index+=1
            try:
                for i in t[tissue_index:]: tissues.append(i)
            except Exception:
                for i in t[1:]: tissues.append(i)
        else:
            try:
                gene = t[0]
                tissue_exp = map(float, t[1:])
                tissue_specific_db[gene]=x,tissue_exp ### Use this to only grab relevant gene expression profiles from the input dataset
            except ZeroDivisionError:
                gene = string.split(t[ens_index],'|')[0] ### Only consider the first listed gene - this gene is the best option based on ExpressionBuilder rankings
                #if 'Pluripotent Stem Cells' in t[marker_in] or 'Heart' in t[marker_in]:
                #if t[marker_in] not in tissues_added: ### Only add the first instance of a gene for that tissue - used more for testing to quickly run the analysis
                tissue_exp = map(float, t[tissue_index:])
                if value_type == 'calls':
                    tissue_exp = produceDetectionCalls(tissue_exp,platform) ### 0 or 1 calls
                tissue_specific_db[gene]=x,tissue_exp ### Use this to only grab relevant gene expression profiles from the input dataset
                tissues_added[t[marker_in]]=[]
            x+=1
    print len(tissue_specific_db), 'genes in the tissue compendium database'

    if correlate_to_tissue_specific == 'yes':
        try: importTissueCorrelations(filename)
        except Exception:
            null=[]
            #print '\nNo tissue-specific correlations file present. Skipping analysis.'; kill
        
def importTissueCorrelations(filename):
    filename = string.replace(filename,'specific','specific_correlations')
    fn=filepath(filename); x=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x==0: x=1 ### Ignore header line
        else:
            uid,symbol,rho,tissue = string.split(data,'\t')
            if float(rho)>rho_threshold: ### Variable used for testing different thresholds internally
                try: tissue_to_gene[tissue].append(uid)
                except Exception: tissue_to_gene[tissue] = [uid]

def simpleUIDImport(filename):
    """Import the UIDs in the gene expression file"""
    uid_db={}
    fn=filepath(filename)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        uid_db[string.split(data,'\t')[0]]=[]
    return uid_db
        
def importGeneExpressionValues(filename,tissue_specific_db,keyed_by,translation_db):
    ### Import gene-level expression raw values           
    fn=filepath(filename); x=0; genes_added={}; gene_expression_db={}
    dataset_name = string.split(filename,delim)[-1][:-4]
    #print 'importing:',dataset_name
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        
        if x==0:
            if '#' not in data:
                for i in t[1:]: sample_headers.append(i)
                x=1
        else:
            gene = t[0]
            #if '-' not in gene and ':E' in gene: print gene;sys.exit()
            if analysis_type == 'AltExon':
                try: ens_gene,exon = string.split(gene,'-')[:2]
                except Exception: exon = gene
                gene = exon
            if keyed_by == 'translation': ### alternative value is 'primaryID'
                """if gene == 'ENSMUSG00000025915-E19.3':
                    for i in translation_db: print [i], len(translation_db); break
                    print gene, [translation_db[gene]];sys.exit()"""
                try: gene = translation_db[gene] ### Ensembl annotations
                except Exception: gene = 'null'
                
            if gene in tissue_specific_db:
                index,tissue_exp=tissue_specific_db[gene]
                try: genes_added[gene]+=1
                except Exception: genes_added[gene]=1
                try: exp_vals = map(float, t[1:])
                except Exception:
                    ### If a non-numeric value in the list
                    exp_vals=[]
                    for i in t[1:]:
                        try: exp_vals.append(float(i))
                        except Exception: exp_vals.append(i)
                if value_type == 'calls': ### Hence, this is a DABG or RNA-Seq expression
                    exp_vals = produceDetectionCalls(exp_vals,targetPlatform) ### 0 or 1 calls
                gene_expression_db[gene] = [index,exp_vals]

    #print len(gene_expression_db), 'matching genes in the dataset and tissue compendium database'
    
    for gene in genes_added:
        if genes_added[gene]>1: del gene_expression_db[gene] ### delete entries that are present in the input set multiple times (not trustworthy)
        else: expession_subset.append(gene_expression_db[gene]) ### These contain the rank order and expression
    #print len(expession_subset);sys.exit()
    expession_subset.sort() ### This order now matches that of 
    gene_expression_db=[]

def produceDetectionCalls(values,Platform):
    # Platform can be the compendium platform (targetPlatform) or analyzed data platform (platform or array_type)
    new=[]
    for value in values:
        if Platform == 'RNASeq':
            if value>1:
                new.append(1) ### expressed
            else:
                new.append(0)
        else:
            if value<cutoff: new.append(1)
            else: new.append(0)
    return new

def importGeneIDTranslations(filename):
    ### Import ExpressionOutput/DATASET file to obtain Ensembl associations (typically for Affymetrix 3' arrays)
    fn=filepath(filename); x=0; translation_db={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if x==0:
            headers = t; x=1; index=0
            for i in headers:
                if 'Ensembl' in i: ens_index = index; break
                index+=1
        else:
            uid = t[0]
            ens_geneids = t[ens_index]
            ens_geneid = string.split(ens_geneids,'|')[0] ### In v.2.0.5, the first ID is the best protein coding candidate
            if len(ens_geneid)>0:
                translation_db[uid] = ens_geneid
    return translation_db

def remoteImportExonIDTranslations(array_type,species,translate_to_genearray,targetplatform):
    global targetPlatform; targetPlatform = targetplatform
    translation_db = importExonIDTranslations(array_type,species,translate_to_genearray)
    return translation_db
    
def importExonIDTranslations(array_type,species,translate_to_genearray):
    gene_translation_db={}; gene_translation_db2={}
    if targetPlatform == 'gene' and translate_to_genearray == 'no':
        ### Get gene array to exon array probeset associations
        gene_translation_db = importExonIDTranslations('gene',species,'yes')
        for geneid in gene_translation_db:
            exonid = gene_translation_db[geneid]
            gene_translation_db2[exonid] = geneid
            #print exonid, geneid
        translation_db = gene_translation_db2
    else:

        filename = 'AltDatabase/'+species+'/'+array_type+'/'+species+'_'+array_type+'-exon_probesets.txt'
        ### Import exon array to target platform translations (built for DomainGraph visualization)
        fn=filepath(filename); x=0; translation_db={}
        print 'Importing the translation file',string.split(fn,delim)[-1][:-4]
        for line in open(fn,'rU').xreadlines():
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            if x==0:  x=1
            else:
                platform_id,exon_id = t
                if targetPlatform == 'gene' and translate_to_genearray == 'no':
                    try:
                        translation_db[platform_id] = gene_translation_db[exon_id] ### return RNA-Seq to gene array probeset ID
                        #print platform_id, exon_id, gene_translation_db[exon_id];sys.exit()
                    except Exception: null=[]
                else:
                    translation_db[platform_id] = exon_id
        del gene_translation_db; del gene_translation_db2
    return translation_db

def analyzeTissueSpecificExpressionPatterns():
    tissue_specific_sorted = []; genes_present={}; tissue_exp_db={}; gene_order_db={}; gene_order=[]
    for (index,vals) in expession_subset: genes_present[index]=[]
    for gene in tissue_specific_db:
        tissue_specific_sorted.append(tissue_specific_db[gene])
        gene_order_db[tissue_specific_db[gene][0]] = gene ### index order (this index was created before filtering)
    tissue_specific_sorted.sort()

    new_index=0    
    for (index,tissue_exp) in tissue_specific_sorted:
        try:
            null=genes_present[index]
            i=0
            gene_order.append([new_index,gene_order_db[index]]); new_index+=1
            for f in tissue_exp:
                ### The order of the tissue specific expression profiles is based on the import gene order
                try: tissue_exp_db[tissues[i]].append(f)
                except Exception: tissue_exp_db[tissues[i]] = [f]
                i+=1
            
        except Exception: null=[] ### Gene is not present in the input dataset

    ### Organize sample expression, with the same gene order as the tissue expression set
    sample_exp_db={}
    for (index,exp_vals) in expession_subset:
        i=0
        for f in exp_vals:
            ### The order of the tissue specific expression profiles is based on the import gene order
            try: sample_exp_db[sample_headers[i]].append(f)
            except Exception: sample_exp_db[sample_headers[i]] = [f]
            i+=1

    if correlate_by_order == 'yes':
        ### Rather than correlate to the absolute expression order, correlate to the order of expression (lowest to highest)
        sample_exp_db = replaceExpressionWithOrder(sample_exp_db)
        tissue_exp_db = replaceExpressionWithOrder(tissue_exp_db)

    global tissue_comparison_scores; tissue_comparison_scores={}
    
    if correlate_to_tissue_specific == 'yes':
        ### Create a gene_index that reflects the current position of each gene
        gene_index={}
        for (i,gene) in gene_order: gene_index[gene] = i
        ### Create a tissue to gene-index from the gene_index 
        tissue_to_index={}
        for tissue in tissue_to_gene:
            for gene in tissue_to_gene[tissue]:
                if gene in gene_index: ### Some are not in both tissue and sample datasets
                    index = gene_index[gene] ### Store by index, since the tissue and expression lists are sorted by index
                    try: tissue_to_index[tissue].append(index)
                    except Exception: tissue_to_index[tissue] = [index]
            tissue_to_index[tissue].sort()
        sample_exp_db,tissue_exp_db = returnTissueSpecificExpressionProfiles(sample_exp_db,tissue_exp_db,tissue_to_index)
        
    PearsonCorrelationAnalysis(sample_exp_db,tissue_exp_db)
    sample_exp_db=[]; tissue_exp_db=[]
    zscore_output_dir,tissue_scores = exportCorrelationResults()
    return zscore_output_dir, tissue_scores

def returnTissueSpecificExpressionProfiles(sample_exp_db,tissue_exp_db,tissue_to_index):
    tissue_exp_db_abreviated={}
    sample_exp_db_abreviated={} ### This db is designed differently than the non-tissue specific (keyed by known tissues)

    ### Build the tissue specific expression profiles    
    for tissue in tissue_exp_db:
        tissue_exp_db_abreviated[tissue] = []
        for index in tissue_to_index[tissue]:
            tissue_exp_db_abreviated[tissue].append(tissue_exp_db[tissue][index]) ### populate with just marker expression profiles

    ### Build the sample specific expression profiles
    for sample in sample_exp_db:
        sample_tissue_exp_db={}
        sample_exp_db[sample]
        for tissue in tissue_to_index:
            sample_tissue_exp_db[tissue] = []
            for index in tissue_to_index[tissue]:
                sample_tissue_exp_db[tissue].append(sample_exp_db[sample][index])
        sample_exp_db_abreviated[sample] = sample_tissue_exp_db
    return sample_exp_db_abreviated, tissue_exp_db_abreviated

def replaceExpressionWithOrder(sample_exp_db):
    for sample in sample_exp_db:
        sample_exp_sorted=[]; i=0
        for exp_val in sample_exp_db[sample]: sample_exp_sorted.append([exp_val,i]); i+=1
        sample_exp_sorted.sort(); sample_exp_resort = []; order = 0
        for (exp_val,i) in sample_exp_sorted: sample_exp_resort.append([i,order]); order+=1
        sample_exp_resort.sort(); sample_exp_sorted=[] ### Order lowest expression to highest
        for (i,o) in sample_exp_resort: sample_exp_sorted.append(o) ### The expression order replaces the expression, in the original order
        sample_exp_db[sample] = sample_exp_sorted ### Replace exp with order
    return sample_exp_db

def PearsonCorrelationAnalysis(sample_exp_db,tissue_exp_db):
    #print "Beginning LineageProfiler analysis"
    k=0
    original_increment = int(len(tissue_exp_db)/15.00); increment = original_increment
    p = 1 ### Default value if not calculated
    for tissue in tissue_exp_db:
        #print k,"of",len(tissue_exp_db),"classifier tissue/cell-types"
        if k == increment: increment+=original_increment; #print '*',
        k+=1
        tissue_expression_list = tissue_exp_db[tissue]
        for sample in sample_exp_db:
            if correlate_to_tissue_specific == 'yes':
                ### Keyed by tissue specific sample profiles
                sample_expression_list = sample_exp_db[sample][tissue] ### dictionary as the value for sample_exp_db[sample]
                #print tissue, sample_expression_list
                #print tissue_expression_list; sys.exit()
            else: sample_expression_list = sample_exp_db[sample]
            try:
                ### p-value is likely useful to report (not supreemly accurate but likely sufficient)
                rho,p = stats.pearsonr(tissue_expression_list,sample_expression_list)
                try: tissue_comparison_scores[tissue].append([rho,p,sample])
                except Exception: tissue_comparison_scores[tissue] = [[rho,p,sample]]
            except Exception:
                ### simple pure python implementation - no scipy required (not as fast though and no p-value)
                #rho = pearson(tissue_expression_list,sample_expression_list)
                None
            #tst = salstat_stats.TwoSampleTests(tissue_expression_list,sample_expression_list)
            #pp,pr = tst.PearsonsCorrelation()
            #sp,sr = tst.SpearmansCorrelation()
            #print tissue, sample
            #if rho>.5: print [rho, pr, sr],[pp,sp];sys.exit()
            #if rho<.5: print [rho, pr, sr],[pp,sp];sys.exit()
    sample_exp_db=[]; tissue_exp_db=[]
    #print 'Correlation analysis finished'
    
def pearson(array1,array2):
    item = 0; sum_a = 0; sum_b = 0; sum_c = 0    
    while item < len(array1):
        a = (array1[item] - avg(array1))*(array2[item] - avg(array2))
        b = math.pow((array1[item] - avg(array1)),2)
        c = math.pow((array2[item] - avg(array2)),2)        
        sum_a = sum_a + a
        sum_b = sum_b + b
        sum_c = sum_c + c
        item = item + 1
    r = sum_a/math.sqrt(sum_b*sum_c)
    return r

def avg(array):
    try: return sum(array)/len(array)
    except Exception: return 0

def adjustPValues():
    """ Can be applied to calculate an FDR p-value on the p-value reported by scipy.
        Currently this method is not employed since the p-values are not sufficiently
        stringent or appropriate for this type of analysis """
        
    import statistics
    all_sample_data={}
    for tissue in tissue_comparison_scores:
        for (r,p,sample) in tissue_comparison_scores[tissue]:
            all_sample_data[sample] = db = {} ### populate this dictionary and create sub-dictionaries
        break
    
    for tissue in tissue_comparison_scores:
        for (r,p,sample) in tissue_comparison_scores[tissue]:
            gs = statistics.GroupStats('','',p)
            all_sample_data[sample][tissue] = gs 
    for sample in all_sample_data:
        statistics.adjustPermuteStats(all_sample_data[sample])
     
    for tissue in tissue_comparison_scores:
        scores = []
        for (r,p,sample) in tissue_comparison_scores[tissue]:
            p = all_sample_data[sample][tissue].AdjP()
            scores.append([r,p,sample])
        tissue_comparison_scores[tissue] = scores

def stdev(array):
    sum_dev = 0
    x_bar = scipy.average(array)
    n = float(len(array))
    for x in array:
        x = float(x)
        sq_deviation = math.pow((x-x_bar),2)
        sum_dev += sq_deviation

    try:
        s_sqr = (1.0/(n-1.0))*sum_dev #s squared is the variance
        s = math.sqrt(s_sqr)
    except ZeroDivisionError:
        s = 'null'
    return s

def replacePearsonPvalueWithZscore():
    all_sample_data={}
    for tissue in tissue_comparison_scores:
        for (r,p,sample) in tissue_comparison_scores[tissue]:
            all_sample_data[sample] = [] ### populate this dictionary and create sub-dictionaries
        break

    for tissue in tissue_comparison_scores:
        for (r,p,sample) in tissue_comparison_scores[tissue]:
            all_sample_data[sample].append(r)

    sample_stats={}
    all_dataset_rho_values=[]
    ### Get average and standard deviation for all sample rho's
    for sample in all_sample_data:
        all_dataset_rho_values+=all_sample_data[sample]
        avg=scipy.average(all_sample_data[sample])
        st_dev=stdev(all_sample_data[sample])
        sample_stats[sample]=avg,st_dev
    
    global_rho_avg = scipy.average(all_dataset_rho_values)
    global_rho_stdev = stdev(all_dataset_rho_values)
    
    ### Replace the p-value for each rho
    for tissue in tissue_comparison_scores:
        scores = []
        for (r,p,sample) in tissue_comparison_scores[tissue]:
            #u,s=sample_stats[sample]
            #z = (r-u)/s
            z = (r-global_rho_avg)/global_rho_stdev ### Instead of doing this for the sample background, do it relative to all analyzed samples
            scores.append([r,z,sample])
        tissue_comparison_scores[tissue] = scores

def exportCorrelationResults():
    corr_output_file = string.replace(exp_output_file,'DATASET','LineageCorrelations')
    corr_output_file = string.replace(corr_output_file,'.txt','-'+coding_type+'.txt')
    if analysis_type == 'AltExon':
        corr_output_file = string.replace(corr_output_file,coding_type,'AltExon')
    filename =  string.split(corr_output_file,delim)[-1][:-4]
    #score_data = exportFile(corr_output_file)
    if use_scipy:
        zscore_output_dir = string.replace(corr_output_file,'.txt','-zscores.txt')
        probability_data = exportFile(zscore_output_dir)
        #adjustPValues()
        replacePearsonPvalueWithZscore()
    ### Make title row
    headers=['Sample_name']
    for tissue in tissue_comparison_scores:
        for (r,z,sample) in tissue_comparison_scores[tissue]: headers.append(sample)
        break
    title_row = string.join(headers,'\t')+'\n'
    #score_data.write(title_row)
    if use_scipy:
        probability_data.write(title_row)
    ### Export correlation data
    tissue_scores = {}; tissue_probabilities={}; tissue_score_list = [] ### store and rank tissues according to max(score)
    for tissue in tissue_comparison_scores:
        scores=[]
        probabilities=[]
        for (r,z,sample) in tissue_comparison_scores[tissue]:
            scores.append(r)
            probabilities.append(z)
        tissue_score_list.append((max(scores),tissue))
        tissue_scores[tissue] = probabilities ### These are actually z-scores
        #tissue_scores[tissue] = string.join(map(str,[tissue]+scores),'\t')+'\n' ### export line
        if use_scipy:
            tissue_probabilities[tissue] = string.join(map(str,[tissue]+probabilities),'\t')+'\n'
        
    tissue_score_list.sort()
    tissue_score_list.reverse()
    #for (score,tissue) in tissue_score_list:
    #score_data.write(tissue_scores[tissue])
    #if use_scipy: probability_data.write(tissue_probabilities[tissue])
    #score_data.close()
    if use_scipy:
        probability_data.close()
    #print filename,'exported...'
    tissue_scores['headers'] = headers
    return zscore_output_dir, tissue_scores

def visualizeLineageZscores(zscore_output_dir,grouped_lineage_zscore_dir,graphic_links):
    import clustering
    ### Perform hierarchical clustering on the LineageProfiler Zscores
    graphic_links = clustering.runHCOnly(zscore_output_dir,graphic_links)   
    return graphic_links
  
if __name__ == '__main__':
    
    ################  Default Variables ################
    species = 'Hs'
    platform = "exon"
    vendor = 'Affymetrix'
    compendium_platform = "exon"
    codingtype = 'protein_coding'
    platform = vendor, platform
    exp_output = None
    geneModels = False
    modelSize = None

    """ This script iterates the LineageProfiler algorithm (correlation based classification method) to identify sample types relative
    two one of two references given one or more gene models. The program '
    """
    #python LineageProfilerIterate.py --i "/Users/nsalomonis/Desktop/dataAnalysis/Sarwal/qPCR/deltaCT/LabMeeting/ExpressionInput/exp.ABI_PediatricSNS.txt" --r "/Users/nsalomonis/Desktop/dataAnalysis/Sarwal/qPCR/deltaCT/LabMeeting/ExpressionOutput/MarkerFinder/MarkerFinder-ABI_PediatricSNS.txt" --m "/Users/nsalomonis/Desktop/dataAnalysis/Sarwal/qPCR/deltaCT/LabMeeting/ExpressionInput/7GeneModels.txt"
    #python LineageProfilerIterate.py --i "/Users/nsalomonis/Desktop/dataAnalysis/Sarwal/qPCR/deltaCT/LabMeeting/ExpressionInput/exp.ABI_PediatricSNS.txt" --r "/Users/nsalomonis/Desktop/dataAnalysis/Sarwal/qPCR/deltaCT/LabMeeting/ExpressionOutput/MarkerFinder/MarkerFinder-ABI_PediatricSNS.txt" --s 4

    ################  Comand-line arguments ################
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Please designate a tab-delimited input expression file in the command-line"
        print 'Example: python LineageProfilerIterate.py --i "/Users/me/qPCR.txt" --r "/Users/me/reference.txt" --m "/Users/me/models.txt"'
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','species=','o=','platform=','codingtype=',
                                                    'compendium_platform=','r=','m=','v=','s='])
        for opt, arg in options:
            if opt == '--i': exp_input=arg
            elif opt == '--o': exp_output=arg
            elif opt == '--platform': platform=arg
            elif opt == '--codingtype': codingtype=arg
            elif opt == '--compendium_platform': compendium_platform=arg
            elif opt == '--r': customMarkers=arg
            elif opt == '--m': geneModels=arg
            elif opt == '--v': vendor=arg
            elif opt == '--s':
                try: modelSize = int(arg)
                except Exception:
                    print 'Please specify a modelSize (e.g., 7-gene model search) as a single integer (e.g., 7)'
                    sys.exit()
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
            
        if exp_output == None: exp_output = exp_input
        
        runLineageProfiler(species,platform,exp_input,exp_output,codingtype,compendium_platform,modelSize=modelSize,customMarkers=customMarkers,geneModels=geneModels)