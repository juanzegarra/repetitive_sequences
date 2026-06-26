#AUTHOR: JUAN ZEGARRA
#DATE: 26/06/2026

#pipeline for assembling COI and per gene sequences by reference by first building a consensus reference sequence, then aligning the built sequence back to the original sequences

from Bio import SeqIO
import subprocess
from multiprocessing import Pool
import re
import os
import glob

#fasta file that can contain mixed genes from NCBI, they must have the gene name in their header
gene_fasta = "genes.fasta"

#ids present in the fasta file
ids = {"COI", "per"}

genes = list(SeqIO.parse(gene_fasta, "fasta"))

#find the interest genes inside the input fasta file, create separate files for each gene 
for id in ids:
    pattern = re.compile(rf"\b{id}\b", re.IGNORECASE)
    with open(f"{id}.fasta", "w") as output_file:
        for gene in genes:
            if pattern.search(gene.description):
                output_file.write(f">{gene.id}\n{gene.seq}\n")

#read files and their identifier are organized in tuples (R1, R2, ID)
genomes = {("/home/rhuanmedeiro/job_soap/Dmeri_cleaned_R1.fastq", "/home/rhuanmedeiro/job_soap/Dmeri_cleaned_R2.fastq", "SER"),
           ("/home/rhuanmedeiro/bowtie_cov_sat/SRR26246248_1_cleaned.fastq", "/home/rhuanmedeiro/bowtie_cov_sat/SRR26246248_2_cleaned.fastq", "ANG"), 
           ("home/paired_DMA/SRR26246249_cleaned_1.fastq", "home/paired_DMA/SRR26246249_cleaned_2.fastq", "DMA")}

#use muscle and emboss to create consensus
def create_consensus(id):
    subprocess.run(["muscle", "-in", f"{id}.fasta", "-out", f"{id}_aln.fasta", "-maxiters", "10"], check=True)
    subprocess.run(["cons", "-sequence", f"{id}_aln.fasta", "-outseq", f"{id}_cons.fasta", "-name", f"{id}_cons"], check=True)
    subprocess.run(["bowtie2-build", f"{id}_cons.fasta", f"{id}_cons"], check=True)

    #assemble using bowtie/samtools for each genome
    for genome in genomes:
        subprocess.run([
            "bowtie2", "-x", f"{id}_cons",
            "-1", genome[0], "-2", genome[1],
            "-S", f"{id}_{genome[2]}.sam",
            "--no-unal", "-p", "10"
        ], check=True)

        with open(f"{id}_{genome[2]}.bam", "wb") as f:
            subprocess.run(
                ["samtools", "view", "-b", f"{id}_{genome[2]}.sam"],
                stdout=f, check=True
            )

        subprocess.run([
            "samtools", "sort",
            "-o", f"{id}_{genome[2]}.sorted.bam",
            f"{id}_{genome[2]}.bam"
        ], check=True)

        subprocess.run([
            "samtools", "consensus",
            "-f", "FASTA",
            "-o", f"{id}_{genome[2]}.fasta",
            f"{id}_{genome[2]}.sorted.bam"
        ], check=True)

        samtools_out = f"{id}_{genome[2]}.fasta"
        temp_file = f"temp_{id}_{genome[2]}.fasta"
        with open(samtools_out) as fin, open(temp_file, "w") as fout:
            fout.writelines(
                f">{genome[2]}\n" if l.startswith(">") else l
                for l in fin
            )
        os.replace(temp_file, samtools_out)


    #clean up temp files
    bam_sam = glob.glob(f"{id}_*am")
    index = glob.glob(f"{id}_*.bt2")
    remove_files = bam_sam + index
    for f in remove_files:
        os.remove(f)

    #group reference sequences and assembled sequences
    with open(f"{id}_combined.fasta", "wb") as out:
        for genome in genomes:
            with open(samtools_out, "rb") as f:
                out.write(f.read())
        with open(f"{id}.fasta", "rb") as f:
            out.write(f.read())
    #align them
    subprocess.run([
        "muscle",
        "-in", f"{id}_combined.fasta",
        "-out", f"{id}_combined_aln.fasta",
        "-maxiters", "10"
    ], check=True)

#multithread
if __name__ == "__main__":
    with Pool(processes=5) as p:
        p.map(create_consensus, ids)
