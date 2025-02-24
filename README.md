# COMP483-Pipeline-Project

## Dependencies
Python 3.10.12  
R 4.4.2  
[kallisto](https://pachterlab.github.io/kallisto/download)  
[sleuth](https://pachterlab.github.io/sleuth/download)  
[bowtie2](https://github.com/BenLangmead/bowtie2)  
[SPAdes](https://github.com/ablab/spades)  
[BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)  

## Required Data
The pipeline takes as input a folder of fastq files and the GenBank accession number for a reference genome  
I used `prefetch` and `fasterq-dump` from [sra-tools](https://github.com/ncbi/sra-tools)  
Test data can be found in `testdata`  
Test data runs with accession number `NC_006273.2`  

## Overview
### Step 1
etc., etc.