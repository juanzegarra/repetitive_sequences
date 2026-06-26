#!/bin/bash
#SBATCH --job-name=bowtie2_sats
#SBATCH --output=out.txt
#SBATCH --error=err.txt
#SBATCH --partition=fila1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=18:00:00
#Juan Zegarra
#Dir structure= main dir with subdirs containing bowtie2 db files, executes search in each 
#subdir. Subdir name must contain db prefix
MAIN_PATH=$(pwd)

run_samtools(){
local genome=$1
for dir in $MAIN_PATH/*/ ; do base_name=$(basename "$dir");
	cd $dir
        samtools view -b -F 4 "${base_name}_${genome}.sam" > "${base_name}_${genome}.bam";
        samtools sort -o "${base_name}_${genome}.sorted.bam" "${base_name}_${genome}.bam";
        samtools index "${base_name}_${genome}.sorted.bam";
        samtools flagstat "${base_name}_${genome}.sorted.bam";
        cd ..;
done
}

for dir in $MAIN_PATH/*/ ; do
        base_name=$(basename "$dir");
        cd $base_name ; 
        srun --mpi=none --ntasks=1 --exact bowtie2 -x $base_name -1 ~/paired_DMA/SRR26246249_cleaned_1.fastq -2 ~/paired_DMA/SRR26246249_cleaned_2.fastq -S "${base_name%/}_DMA.sam" --no-unal -p 15 &
        srun --mpi=none --ntasks=1 --exact bowtie2 -x $base_name -1 ~/bowtie_cov_sat/SRR26246248_1_cleaned.fastq -2 ~/bowtie_cov_sat/SRR26246248_2_cleaned.fastq -S "${base_name%/}_ANG.sam" --no-unal -p 15 &
        srun --mpi=none --ntasks=1 --exact bowtie2 -x coi_dmer_cons.fasta -1 ~/job_soap/Dmeri_cleaned_R1.fastq -2 ~/job_soap/Dmeri_cleaned_R2.fastq -S "${base_name%/}_SER.sam" --no-unal -p 15&
        wait
        cd ..;
done

genomes={"SER" "ANG" "DMA"}

for i in "${genomes[@]}" do; 
	run_samtools "${i}"
