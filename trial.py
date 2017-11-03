import pysam
import string,os,sys,copy
import time
import getopt

def findGeneVariants(species,symbols,bam_dir,variants=None):
    global insertion_db
    insertion_db={}
    
    if len(symbols)>0:
        ### Search for genes not for coordinates
        search_locations = geneCoordinates(species,symbols)
    else:
        ### Search for coordinates and not genes
        search_locations = variantCoordinates(variants)
        
    print search_locations
    ### Discover the variants
    #variant_db = findVariants(bam_dir,search_locations)
    
    #variant_filtered_db={}
    #for var in variant_db:
        #print var, variant_db[var]
        #if variant_db[var]>3:
            #print var,variant_db[var]
           # variant_filtered_db[var] = variant_db[var]

def geneCoordinates(species,symbols):
    genes=[]
    import EnsemblImport
    ensembl_annotation_db = EnsemblImport.reimportEnsemblAnnotations(species,symbolKey=True)
    for symbol in symbols:
        if symbol in ensembl_annotation_db:
            ens_geneid = ensembl_annotation_db[symbol]
            #ensembl_annotation_db[ensembl_gene_id] = symbol,description,mRNA_processing
            genes.append((ens_geneid,symbol))
        else:
            print symbol, 'not found'


if __name__ == "__main__":
    #bam_dir = "H9.102.2.6.bam"
    #reference_dir = 'H9.102.2.6__exon.bed'
    ################  Comand-line arguments ################
    symbols=[]
    variantFile = None
    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print "Warning! Please designate a BAM file as input in the command-line"
        print "Example: python BAMtoExonBED.py --i /Users/me/sample1.bam --r /Users/me/Hs_exon-cancer_hg19.bed"
        sys.exit()
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','species=','g=','v='])
        #both options and remiander have been assigned the values(same)
        for opt, arg in options:
            #http://www.tutorialspoint.com/python/python_command_line_arguments.htm
            if opt == '--i': bam_dir=arg ### A single BAM file location (full path)
            elif opt == '--species': species=arg 
            elif opt == '--g': symbols.append(arg)
            elif opt == '--v': variantFile = arg
            else:
                print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
    findGeneVariants(species,symbols,bam_dir,variants=variantFile)
