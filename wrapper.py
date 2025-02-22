import os
import glob
import subprocess
import re
import pandas as pd

    ### EFETCH FOR TARGET ORGANISM SOMEHOW MAYBE ###

def download_cds(acc):
    os.mkdir("temp") # store ncbi datasets temporarily
    os.chdir("temp")
    os.system(f"datasets download virus genome accession {acc} --include cds") # get all protein sequences for herpesvirus
    os.system("unzip ncbi_dataset.zip")
    os.system(f"mv ncbi_dataset/data/cds.fna ../herpes.fasta") # get fasta file
    os.chdir("..")
    os.system("rm -r temp") # get rid of everything else
    return

def count_cds(log):
    with open("herpes.fasta") as f: # count number of CDS
        cds = len(re.findall(">", f.read())) # use regedx matching to find individual cds
    with open(log, "w") as log:
        log.write("The HCMV genome (NC_006273.2) has {str(cds)} CDS\n\n")
    return

def run_kallisto(ind):
    os.system(f"kallisto index -i {ind} herpes.fasta") # build kallisto index from protein sequences
    os.mkdir("results")
    paths = sorted(glob.glob("../data/fastq_files/*")) # get all fastq files
    for i in range(0, len(paths), 2): # get each pair of files
        pair1 = paths[i]
        pair2 = paths[i+1]
        tag = pair1.split("/")[-1].split("_")[0] # get the srr acc for the output dir
        os.mkdir(f"results/{tag}") # make output dir to store results for each pair
        os.system(f"kallisto quant -i {ind} -o results/{tag} -b 10 -t 2 {pair1} {pair2}") # run kallisto on each pair
    return

def write_tpm_stats(log):
    header = ["sample", "condition", "min_tpm", "med_tpm", "mean_tpm", "max_tpm"]
    samples = ["Donor 1", "Donor 1", "Donor 3", "Donor 3"]
    conditions = ["2dpi", "6dpi", "2dpi", "6dpi"]
    with open(log, "a") as log:
        log.write("\t".join(header) + "\n")
        results = sorted(glob.glob("results/*"))
        for i in range(len(results)):
            r = results[i]
            tpm = pd.read_csv(f"{r}/abundance.tsv", sep="\t")["tpm"]
            row = [samples[i], conditions[i], str(tpm.min()), str(tpm.median()), str(tpm.mean()), str(tpm.max())]
            log.write("\t".join(row) + "\n")
        log.write("\n")
    return

def run_sleuth():
    results = sorted(glob.glob("results/*"))
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

def main():
    os.chdir("PipelineProject_Ben_Moginot")

    log = "PipelineProject.log"

    acc = "NC_006273.2"
    download_cds(acc)
    count_cds(log)

    index = "herpes.idx"
    run_kallisto(index)

    write_tpm_stats(log)

    run_sleuth()
    write_diff_exp(log)

if __name__ == "__main__":
    main()