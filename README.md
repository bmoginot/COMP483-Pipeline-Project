# COMP483-Pipeline-Project

## Dependencies
Python 3.10.12  
R 4.4.2  
[BioPython](https://biopython.org/wiki/Download)  
[pandas](https://pandas.pydata.org/docs/getting_started/install.html)  
[sleuth](https://pachterlab.github.io/sleuth/download)    
[kallisto](https://pachterlab.github.io/kallisto/download)  
[bowtie2](https://github.com/BenLangmead/bowtie2)  
[SPAdes](https://github.com/ablab/spades)  
[BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)  

## Required Data
This pipeline takes a folder of paired-end fastq files as input  
I used `prefetch` and `fasterq-dump` from [sra-tools](https://github.com/ncbi/sra-tools) to retrieve the test data  
Test data can be found in `testdata`   

## How to Run
After cloning the repo, move into it via
```
cd COMP483-Pipeline-Project
```
Run the pipeline with
```
python wrapper.py -i input_data_dir
```
For example, to run the test data, use
```
cd COMP483-Pipeline-Project
python wrapper.py -i testdata
```
