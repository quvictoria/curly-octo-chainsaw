# Frequency analysis for PE sequencing
#### Author: Victoria Qu

## Overview
The pipeline outlined in `deletion_frequency.py` performs a simple deletion frequency calculation for paired end sequencing. The intended purpose is to quantify deletions after editing via a Cas9 nuclease. Ths pipeline utilizes the `BBmap` ([ref](https://anaconda.org/bioconda/bbmap)) tool suite for read merging and alignment and the `pysam` ([ref](https://pysam.readthedocs.io/en/latest/installation.html)) package for SAM/BAM file manipulation. An included `Dockerfile` also lists necessary dependencies for the script. 

## Methods
Inputs: 
- `in_fastqR1`: R1 input FASTQ file (gzipped)
- `n_fastqR2`: R2 input FASTQ file (gzipped)
- `output_dir`: name of desired output directory
- `sample_name`: name of sample, used to generate intermediate and output files
- `ref_fasta`: Reference genome (FA file) 
- `chr`: specified region chromosome (eg `CHR1`)
- `region_start`: starting coordinate of specified region, specified as an integer
- `region_end`: ending coordinate of specified region, specified as an integer
  
Steps:
1. Merging reads (`merge_reads()`): Takes `in_fastqR1` and `in_fastqR2` and merges reads together using `bbmerge`. This outputs a merged FASTQ file that will be fed into the following step. 
2. Count merged reads (`count_merge_reads()`): Counts # of reads in the resulting merged FASTQ file. 
3. Aligning reads (`align_reads()`): Takes merged FASTQ files and aligns to the specified genome reference FA file using `bbmap`. Creates intermediate SAM file, converts SAM to BAM format, sorts and indexes BAM file. 
4. Count aligned reads (`count_aligned_reads()`): Counts # of reads in the intermediate SAM file. 
5. Count reads of interest (`count_reads_in_region()`):  Obtains counts of aligned reads within target region, aligned reads within target region that contain deletions (`contains_deletion()`), and lengths of deletions within target regions. To see if a read contains a deletion.
6. Calculate deletion frequency: (# of reads containing deletions / # of reads within region) *100. 
7. Plot and save histogram of deletion lengths for sample. 

## Future Work
The approach presented here is rather simple, as it only identifies deletions within reads but ignores other editing results such as insertions or substitutions and neglects storage of deletion sequences. Expanding on the currently existing functions in the script, I could add additional logic to capture insertions and substitutions. This would provide a more holistic view of the editing implications for downstream analysis.

For an even more complex approach with inspiration from Crispresso2 from the Pinello lab ([ref](https://github.com/pinellolab/CRISPResso2)), I would generate an allele frequency table of all aligned reads. To do so, we can: 
- Add in additional parameters for amplicon sequence instead of target region coordinates.
- For each read, convert `CIGAR` string into a tuple of uniform length (padding if necessary) that contains the reference sequence and the aligned read sequence.
- From there, I could create resulting nucleotide frequency tables and modification frequency tables that will allow me to customize what fraction of the specified region is included in quantification. That would be more customizable with other Cas nuclease types with different editing signatures. 
