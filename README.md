# COMP483-Pipeline-Project

## Download Data
I used `prefetch` and `fasterq-dump` from [sra-tools](https://github.com/ncbi/sra-tools)

## How it works
build database of transcriptome for target organism (need to take this as input)
output how many coding regions are in target organism > log
run each pair of reads through kallisto and get transcripts per million (TPM)
get stats on TPM > log