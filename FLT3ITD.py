#!/usr/bin/env python

import pysam
import string,os,sys,copy
import time
import getopt


def findVariants(bam_dir,search_locations,multi=False):

  for fi in os.listdir(bam_dir):
    
   if '.bam' in fi and ".bai" not in fi:
    
    bamfile = pysam.Samfile(os.path.join(bam_dir, fi), "rb" )
    
    #output_bed_rows=0
    #https://www.biostars.org/p/76119/
    #variant_db={}
    reference_rows=0

    o = open (string.replace(os.path.join(bam_dir, fi),'.bam','__flt3reads.txt'),"w")
    for (chr,strand,start,stop,symbol) in search_locations: ### read each line one-at-a-time rather than loading all in memory
            #read_count=0
            #reference_rows+=1
            stop=int(stop)+100 ### buffer for single variants
            start=int(start)-100 ### buffer for single variants
            for alignedread in bamfile.fetch(chr, int(start),int(stop)):
                #md = alignedread.opt('MD')
               
                o.write(str(alignedread))
                o.write(str(alignedread.pos)+'\t')
                o.write(alignedread.seq+'\t')
                o.write('\n')
                try:
                    mate = bamfile.mate(alignedread) #https://groups.google.com/forum/#!topic/pysam-user-group/9HM6nx_f2CI
                    o.write(str(mate))
                except Exception:
                     o.write("mate is unmapped")
                o.write('\n')
               
                #omd = md
                
#                codes = map(lambda x: x[0],alignedread.cigar)
#                cigarstring = alignedread.cigarstring
#                #print symbol,cigarstring
#                if 1 in codes and alignedread.pos:
#                    ### Thus an insertion is present
#                    cigarstring = alignedread.cigarstring
#                    chr = bamfile.getrname(alignedread.rname)
#                    pos = alignedread.pos
#                    def getInsertions(cigarList,X):
#                        cummulative=0
#                        coordinates=[]
#                        for (code,seqlen) in cigarList:
#                            if code == 0 or code == 3:
#                                cummulative+=seqlen
#                            if code == 1:
#                                coordinates.append(X+cummulative)
#                        return coordinates
#                    coordinates = getInsertions(alignedread.cigar,pos)
#                    """
#                    print pos
#                    print coordinates
#                    print alignedread.seq
#                    print codes
#                    print alignedread.cigar
#                    print cigarstring
#                    
#                    print md;sys.exit()
#                    """
#                    for pos in coordinates:
#                        try: variant_db[chr,pos,symbol]+=1
#                        except Exception: variant_db[chr,pos,symbol] = 1
#                        insertion_db[chr,pos]=[]
#                    continue
#                    
#                try:
#                    int(md) ### If an integer, no mismatches or deletions present
#                    continue
#                except Exception:
#                    #print chr, int(start),int(stop)
#                    #print alignedread.get_reference_sequence()
#                    #print alignedread.seq
#                    md = string.replace(md,'C','A')
#                    md = string.replace(md,'G','A')
#                    md = string.replace(md,'T','A')
#                    md = string.split(md,'A')
#                    pos = alignedread.pos
#                    chr = bamfile.getrname(alignedread.rname)
#                    #if omd == '34^GA16': print md, pos
#                    for i in md[:-1]:
#                        try:
#                            pos+=int(i)+1
#                        except Exception:
#                            if i == '':
#                                pos+=+1
#                            elif '^' in i: ### position is equal to the last position
#                                pos+=int(string.split(i,'^')[0])+1
#                                #pass
#                        #if 'CGGATCC' in alignedread.seq: print string.split(alignedread.seq,'CGGATCC')[1],[pos]
#                        try: variant_db[chr,pos,symbol]+=1
#                        except Exception: variant_db[chr,pos,symbol] = 1
#                        
#                    #codes = map(lambda x: x[0],alignedread.cigar)
#            output_bed_rows+=1
#    o.close()
#    bamfile.close()
#    if multi==False:
#        print time.time()-start_time, 'seconds to assign reads for %d entries from %d reference entries' % (output_bed_rows,reference_rows)
#    #print variant_db;sys.exit()
#    return variant_db
#
#def pileupAnalysis(bam_dir,search_locations,multi=False,findallflag=False):
#    start_time = time.time()    
#    bamfile = pysam.Samfile(bam_dir, "rb" )
#    reference_rows=0
#    output_bed_rows=0
#    #https://www.biostars.org/p/76119/
#    variant_db={}
#    o = open (string.replace(bam_dir,'.bam','__variant.txt'),"w")
#    entries = ['chr','position','rare-allele frq','type','depth','gene']
#    o.write(string.join(entries,'\t')+'\n')
#    #print 'Analyzing',len(search_locations),'variants'
#    for (chr,pos,symbol) in search_locations: ### read each line one-at-a-time rather than loading all in memory
#        pos = int(pos)
#        read_count=0
#        reference_rows+=1
#        nucleotide_frequency={}
#        for pileupcolumn in bamfile.pileup(chr,pos,pos+1):
#            # Skip columns outside desired range
#            #print pos, pileupcolumn.pos, pileupcolumn.cigarstring, pileupcolumn.alignment.pos
#            if pileupcolumn.pos == (pos-1):
#                for pileupread in pileupcolumn.pileups:
#                    try: nt = pileupread.alignment.query_sequence[pileupread.query_position]
#                    except Exception,e:
#                        if 'D' in pileupread.alignment.cigarstring:
#                            nt = 'del'
#                        else:
#                            nt = 'ins'
#                    try: nucleotide_frequency[nt]+=1
#                    except Exception: nucleotide_frequency[nt]=1
#        nt_freq_list=[]
#        for nt in nucleotide_frequency:
#            nt_freq_list.append(nucleotide_frequency[nt])
#        s = sum(nt_freq_list)
#        nt_freq_list.sort()
#        
#        try:
#            frq = float(search_locations[chr,pos,symbol])/s ### This fixes that (number of insertions from before)
#        except Exception: frq = '1.000000'; print symbol, pos, nucleotide_frequency, search_locations[chr,pos,symbol]
#        if (chr,pos) in insertion_db:
#            #print 'insertion', chr, pos
#            call = 'insertion'
#            ### For insertions if the inserted base matches the reference base, incorrect freq will be reported
#        elif 'del' in nucleotide_frequency:
#            #frq = float(nt_freq_list[-2])/s
#            call = 'del'
#        else:
#            #frq = float(nt_freq_list[-2])/s
#            call = 'mismatch'
#        if len(nt_freq_list)>1 or call == 'insertion':
#            if frq>0.01 and findallflag==False:
#                frq = str(frq)[:4]
#                entries = [chr,str(pos),str(frq),call,str(s),symbol]
#                o.write(string.join(entries,'\t')+'\n')
#                output_bed_rows+=1
#            else:
#                frq = str(frq)
#                entries = [chr,str(pos),str(frq),call,str(s),symbol]
#                o.write(string.join(entries,'\t')+'\n')
#                output_bed_rows+=1
#    
#    o.close()
#    bamfile.close()
#    if multi==False:
#        print time.time()-start_time, 'seconds to assign reads for %d entries from %d reference entries' % (output_bed_rows,reference_rows)
search_locations=[]
search_locations.append(['chr13','-','28608000','28608600','FLT3'])
bam_dir="/Users/meenakshi/Desktop/mounter2/bams_TCGA_AML"
findVariants(bam_dir,search_locations,multi=False)