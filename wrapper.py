import os
import glob
import subprocess
import re
import pandas as pd
from Bio import SeqIO
import argparse
import sys

def check_arg(args=None):
    parser = argparse.ArgumentParser(description="run pipeline")
    parser.add_argument("-i", "--input",
    help="input directory containing paired-end fastq files",
    required=True)
    return parser.parse_args(args)

def download_cds(acc):
    os.mkdir("tmp") # store ncbi datasets temporarily
    os.chdir("tmp")
    os.system(f"datasets download virus genome accession {acc} --include cds,genome") # get all protein sequences for herpesvirus
    os.system("unzip ncbi_dataset.zip")
    os.system(f"mv ncbi_dataset/data/cds.fna ../HCMV_cds.fq") # get cds
    os.system(f"mv ncbi_dataset/data/genomic.fna ../HCMV_ref.fq") # get reference
    os.chdir("..")
    os.system("rm -r tmp") # get rid of everything else
    return

def count_cds(log):
    with open("HCMV_cds.fq") as f: # count number of CDS
        cds = len(re.findall(">", f.read())) # use regedx matching to find individual cds
    with open(log, "w") as l:
        l.write(f"The HCMV genome (NC_006273.2) has {str(cds)} CDS\n\n")
    return

def run_kallisto(ind, data):
    index = f"{ind}.idx"
    os.system(f"kallisto index -i {index} HCMV_cds.fq") # build kallisto index from protein sequences
    os.mkdir("kallisto_results")
    paths = sorted(glob.glob(f"{data}/*")) # get all fastq files
    for i in range(0, len(paths), 2): # get each pair of files
        pair1 = paths[i]
        pair2 = paths[i+1]
        tag = pair1.split("/")[-1].split("_")[0] # get the srr acc for the output dir
        os.mkdir(f"kallisto_results/{tag}") # make output dir to store results for each pair
        os.system(f"kallisto quant -i {index} -o kallisto_results/{tag} -b 10 -t 2 {pair1} {pair2}") # run kallisto on each pair
    return

def write_tpm_stats(log):
    header = ["sample", "condition", "min_tpm", "med_tpm", "mean_tpm", "max_tpm"]
    samples = ["Donor 1", "Donor 1", "Donor 3", "Donor 3"]
    conditions = ["2dpi", "6dpi", "2dpi", "6dpi"]
    with open(log, "a") as l:
        l.write("\t".join(header) + "\n")
        results = sorted(glob.glob("kallisto_results/*"))
        for i in range(len(results)):
            r = results[i]
            tpm = pd.read_csv(f"{r}/abundance.tsv", sep="\t")["tpm"]
            row = [samples[i], conditions[i], str(tpm.min()), str(tpm.median()), str(tpm.mean()), str(tpm.max())]
            l.write("\t".join(row) + "\n")
        l.write("\n")
    return

def run_sleuth():
    results = sorted(glob.glob("kallisto_results/*"))
    conditions = ["2dpi", "6dpi", "2dpi", "6dpi"]
    header = ["sample", "condition", "path"]
    with open("sleuth_input.txt", "w") as out:
        out.write("\t".join(header) + "\n")
        for i in range(len(results)):
            path = results[i]
            sample = path.split("/")[-1]
            out.write("\t".join([sample, conditions[i], path]) + "\n")
    os.system("Rscript ../run_sleuth.R")

def write_diff_exp(log):
    df = pd.read_csv("sleuth_results.txt", sep=" ")
    out = pd.DataFrame(data={"target_id": df["target_id"],
    "test_stat": df["test_stat"],
    "pval": df["pval"],
    "qval": df["qval"]})
    out.to_csv(log, mode="a", sep="\t", index=False)
    with open(log, "a") as l:
        l.write("\n")

def run_bowtie(ind, data):
    os.mkdir("bowtie_results")
    os.chdir("bowtie_results")
    os.system(f"bowtie2-build ../HCMV_ref.fq {ind}") # build bowtie reference
    paths = sorted(glob.glob(f"{data}/*"))
    for i in range(0, len(paths), 2): # get each pair of files
        pair1 = paths[i]
        pair2 = paths[i+1]
        tag = pair1.split("/")[-1].split("_")[0] # get the srr acc for the output sam
        os.system(f"bowtie2 --quiet -x {ind} -1 {pair1} -2 {pair2} -S {tag}map.sam --al-conc {tag}_mapped.fq") # run bowtie on each pair
    os.chdir("..")

def compare_read_counts(log, data):
    donors = ["Donor 1", "Donor 1", "Donor 3", "Donor 3"]
    conditions = ["2dpi", "6dpi", "2dpi", "6dpi"]
    unfilt = sorted(glob.glob(f"{data}/*_1.fq"))
    mapped = sorted(glob.glob("bowtie_results/*.1.fq"))
    for i in range(len(mapped)):
        with open(unfilt[i]) as f:
            before = int(len(f.readlines()) / 4)
        with open(mapped[i]) as g:
            after = int(len(g.readlines()) / 4)
        with open(log, "a") as l:
            l.write(f"{donors[i]} ({conditions[i]}) had {before} read pairs before Bowtie2 filtering and {after} read pairs after.\n")
    with open(log, "a") as l:
        l.write("\n")

def run_spades(log):
    os.mkdir("spades_results")
    os.chdir("spades_results")
    log = f"../{log}"
    results = sorted(glob.glob("../bowtie_results/*.fq"))
    for i in range(0, len(results), 4): # get each pair of files
        reads = [results[i], results[i+1], results[i+2], results[i+3]]
        tag = reads[0].split("/")[-1].split("_")[0][:-1] + "X"
        command = f"spades.py -k 77 -t 2 --only-assembler --pe-1 1 {reads[0]} --pe-2 1 {reads[1]} --pe-1 2 {reads[2]} --pe-2 2 {reads[3]} -o {tag}"
        with open(log, "a") as l:
            l.write(f"{command}\n")
        os.system(command)
    with open(log, "a") as l:
        l.write("\n")
    os.chdir("..")

def run_blast():
    os.mkdir("blast_results")
    os.chdir("blast_results")

    db = "betaherpesvirinae"
    ref = "all_herpes.fa"

    os.mkdir("tmp") # store ncbi datasets temporarily
    os.chdir("tmp")
    os.system(f"datasets download virus genome taxon {db} --refseq --include genome") # get all protein sequences for herpesvirus
    os.system("unzip ncbi_dataset.zip")
    os.system(f"mv ncbi_dataset/data/genomic.fna ../{ref}") # get cds
    os.chdir("..")
    os.system("rm -r tmp") # get rid of everything else

    os.system(f"makeblastdb -in {ref} -out {db} -title {db} -dbtype nucl")

    results = sorted(glob.glob("../spades_results/*"))

    for res in results:
        for i, record in enumerate(SeqIO.parse(f"{res}/contigs.fasta", "fasta")): # get first contig, which is the longest
            longest = f">{str(record.id)}\n{str(record.seq)}"
            if i == 0:
                break
        tag = res.split("/")[-1]
        with open(f"{tag}_longest.fa", "w") as f:
            f.write(longest)

    contigs = sorted(glob.glob("SRR*"))
    for infile in contigs:
        tag = infile.split("_")[0]
        os.system(f"blastn -query {infile} -db {db} -out {tag}_out.csv -outfmt '10 sacc pident length qstart qend sstart send bitscore evalue stitle'")

def parse_blast(log):
    log = f"../{log}"

    results = sorted(glob.glob("*out.csv"))

    samples = ["Donor 1", "Donor 3"]

    for i in range(len(results)):
        rows = []
        with open(results[i]) as f:
            hits = [line.strip() for line in f.readlines()]
            for h in hits:
                rows.append(h.split(","))

        with open(log, "a") as l:
            header = ["sacc", "pident", "length", "qstart", "qend", "sstart", "send", "bitscore", "evalue", "stitle"]
            l.write(samples[i] + "\n")
            l.write("\t".join(header) + "\n")
            for r in rows:
                l.write("\t".join(r) + "\n")
            l.write("\n")


def main():
    args = check_arg(sys.argv[1:])
    cwd = os.getcwd()
    data = os.path.join(cwd, args.input)

    outdir = "PipelineProject_Ben_Moginot"
    if os.path.isdir(outdir):
        os.system(f"rm -r {outdir}")
    os.mkdir(outdir)
    os.chdir(outdir)

    log = "PipelineProject.log"

    acc = "NC_006273.2"
    ind = "HCMV"

    download_cds(acc)
    count_cds(log)

    run_kallisto(ind, data)
    write_tpm_stats(log)

    run_sleuth()
    write_diff_exp(log)

    run_bowtie(ind, data)
    compare_read_counts(log, data)

    run_spades(log)

    run_blast()
    parse_blast(log)

if __name__ == "__main__":
    main()