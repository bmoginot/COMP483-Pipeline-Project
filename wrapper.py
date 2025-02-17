import os
import glob

os.chdir("./PipelineProject_Ben_Moginot")

log = open("PipelineProject.log", "w")

acc = "NC_006273.2"
os.mkdir("temp") # store ncbi datasets temporarily
os.chdir("temp")
os.system(f"datasets download gene accession {acc}") # get all protein sequences for herpesvirus
os.system("unzip ncbi_dataset.zip")
os.system("mv ncbi_dataset/data/protein.faa ../herpes.fasta") # get fasta file
os.chdir("..")
os.system("rm -r temp") # get rid of everything else

index = "herpes.idx"
os.system(f"kallisto index -i {index} herpes.fasta") # build kallisto index from protein sequences

os.mkdir("results")
paths = sorted(glob.glob("../data/fastq_files/*")) # get all fastq files
for i in range(0, len(paths), 2): # get each pair of files
    pair1 = paths[i]
    pair2 = paths[i+1]
    print(pair1, pair2)
    tag = pair1.split("/")[-1].split("_")[0] # get the srr acc for the output dir
    print(tag)
    os.mkdir(f"results/{tag}") # make output dir to store results for each pair
    os.system(f"kallisto quant -i {index} -o results/{tag} -b 10 -t 2 {pair1} {pair2}") # run kallisto on each pair

log.close()