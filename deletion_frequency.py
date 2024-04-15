""" 
Author: Victoria Qu
Purpose: This script does the following steps for deletion quantification of paired end amplicon sequencing. 
It performs the following steps:
    - Merging of reads using BBMerge (bbmerge.sh)
    - Aligns to specified reference genome using BBMap (bbmap.sh)
    - Focuses on a specific target site where deletions are expected and obtains number of reads that contain deletions 
    - Generates deletion size histogram 
"""
import argparse 
import os 
import sys
import re
from subprocess import Popen, PIPE 
import pysam 
import gzip
import pandas as pd 
import numpy as np 
from Bio import SeqIO
import matplotlib.pyplot as plt 

#############################################################functions that we will use ########################################################################

def merge_reads(in_fastqR1, in_fastqR2, merged_fastq):
    """
    Use bbmerge.sh to merge R1 and R2 fastqs. 

    Params: 
    - in_fastqR1 = R1 fastq file path
    - in_fastqR2 = R2 fastq file path 
    - merged_fastq = resulting merged fastq output path 
    """
    #process = Popen(["bash", "bbmap/bbmerge.sh", "in=" + in_fastqR1, "in2="+in_fastqR2, "out=" + merged_fastq, "maxloose", "qtrim=t", "trimq=11"],
    process = Popen(["bbmerge.sh", "in=" + in_fastqR1, "in2="+in_fastqR2, "out=" + merged_fastq, "maxloose", "qtrim=t", "trimq=11"],
    stdout = PIPE, stderr = PIPE)

    stdout, stderr = process.communicate()
    print("Fastqs Merged: ", merged_fastq)

def count_merge_reads(merged_fastq):
    #count merged reads 
    merge_count = 0

    with gzip.open(merged_fastq, 'rt') as file:
        for _ in SeqIO.parse(file, "fastq"):
            merge_count+=1 
    
    return merge_count

def align_reads(input_merged_fastq, ref_fasta, out_sam, out_bam, sorted_bam):
    """
    Use bbmap.sh to align input merged fastq to refernece FASTA. Use pysam to convert SAM TO BAM.  

    Params: 
    - input_merged_fastq = merged fastq fle path
    - ref_fastq = reference file fasta file path 
    - out_sam = aligned SAM file
    - out_bam = converted BAM (from SAM) file 
    - sorted_bam = sorted/indexed output BAM
    """
    #process = Popen(["bash", "bbmap/bbmap.sh", "ref=" + ref_fasta, "in=" + input_merged_fastq, "out=" + out_sam, "nodisk"],
    process = Popen(["bbmap.sh", "ref=" + ref_fasta, "in=" + input_merged_fastq, "out=" + out_sam, "nodisk"],
    stdout = PIPE, stderr = PIPE)

    stdout, stderr = process.communicate()
    print("SAM created: ", out_sam)

    #sam to bam 
    pysam.view("-bS", "-o", out_bam, out_sam, catch_stdout = False)
    print("SAM to BAM:", out_bam)
    pysam.sort("-o", sorted_bam, out_bam)
    print("Sorted BAM:", sorted_bam)
    pysam.index(sorted_bam)
    print("Indexed BAM:", sorted_bam)

def count_aligned_reads(out_sam):
    #count aligned reads from SAM file
    align_count = 0

    with pysam.AlignmentFile(out_sam, "r") as file:
        for read in file:
            if not read.is_unmapped:
                align_count +=1
    return align_count

def contains_deletion(read, region_start, region_end):
    """
    From read, if CIGAR string contains >=1 deletion within specific region, it contains deletion.

    - read
    - region_start
    - region_end 

    """
    reference_position = read.reference_start #initialize with 0 based leftmost pos
    CIGAR = read.cigarstring
    del_found = False 
    num_dels_in_read = 0
    deletion_locs = []

    for match in re.finditer(r'([MIDNSHP=X])(\d+)', CIGAR):
        #instead of cigar tuples, use re to get matches to get operation and length
        operation = int(match.group(1))
        length = int(match.group(2))

        if operation == 'I': 
            pass #ignore since not reflected in reference positon
        elif operation in['=', 'X','M']:
            #if operation is a mismatch, sequence match, or match skip over it 
            reference_position += length 
        elif operation == 'D':
            #if operation is a deletion then keep track of how long deletion length is 
            del_start = reference_position 
            del_end = reference_position + length 

            if ((del_end > region_start) and (del_start < region_end)):
                #check if deletion is within bounds of specified regon
                del_found = True 
                num_dels_in_read += length

                ##get deletion information here using del_start and del_end from the read sequence 

            reference_position += length 
    
    return del_found, num_dels_in_read


def count_reads_in_region(sorted_bam, chr, region_start, region_end):
    """
    For each read within the BAM file, get total reads that fall within specified region, 
    also get the total reads that contain deletions (ignoring insertions and substitutons)
    """
    reads_in_region = 0 #number of reads that fall within specified region
    deletions_in_region = 0 #number of reads that contain >=1 deletion within specified region
    del_lengths = []

    with pysam.AlignmentFile(sorted_bam, 'rb') as bam_file:
        for read in bam_file.fetch(chr, region_start, region_end):
            if not read.is_unmapped:
                reads_in_region += 1    
                del_found, num_dels_in_read = contains_deletion(read, region_start, region_end)
                if del_found:
                    deletions_in_region +=1
                    del_lengths.append(num_dels_in_read)
                

    return reads_in_region, deletions_in_region, del_lengths

#################################################### CL args #################################################################
parser = argparse.ArgumentParser(description= "Calculate deletion frequency from PE fastqs")
#### add arugments to run in command line for output_dir, fastq R1, fastaR2, ref fasta, sample name, specified region coordinates
#parser.add_argument() 

if __name__ == '__main__':
    output_dir = 'test_output' #add to parser argument 

    #Merge R1 and R2 input fastqs using merge_reads 
    in_fastqR1 = 'TestSample_L001_R1_001.fastq.gz'  #add to parser argument 
    in_fastqR2 = 'TestSample_L001_R2_001.fastq.gz'  #add to parser argument 

    sample_name = 'TestSample'  #add to parser argument 
    merged_fastq = sample_name + 'merged.fastq.gz'  


    merge_reads(in_fastqR1, in_fastqR2, merged_fastq)
    merge_count = count_merge_reads(merged_fastq) #get # of merged reads 

    #Align merged fastq to reference, convert SAM to BAM 
    ref_fasta = 'hs38DH.fa' #add to parser argument 
    out_sam = sample_name + '.sam'
    out_bam = sample_name + '.unsorted.bam'
    sorted_bam = sample_name + '.sorted.bam'

    align_reads(merged_fastq, ref_fasta, out_sam, out_bam, sorted_bam)
    align_count = count_aligned_reads(out_sam) #get # of aligned reads 

    #filter reads to desired region for EMX1 
    chr = 2 #add to parse argument 
    region_start = 72916260 #add to parser argument 
    region_end = 72936071 #add to parser argument 

    #count target reads within region and for reads with deletions within same region 
    reads_within_region, deletions_in_region, del_lengths = count_reads_in_region(sorted_bam, chr, region_start, region_end)

    #get a deletion size histogram to see how the dstibuion of deletion sizes look likes for the BAM alignment 
    deletion_length_df = pd.DataFrame({'deletion_lengths': del_lengths})
    deletion_length_df.hist(bins=3)
    plt.savefig(f'{output_dir}/deletion_hist_for_{sample_name}.png')

    #deletion percentage 
    del_frequency = (deletions_in_region/reads_within_region)*100
    print(f"deletion frequency for {sample_name}:",del_frequency) 
    
    







