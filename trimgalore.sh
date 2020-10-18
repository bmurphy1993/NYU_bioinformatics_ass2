#!/bin/bash
#SBATCH --job-name=NGS1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bem6982@nyu.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/scratch/bem6982/bioinformatics/assn_02_ngs/NGS1_%j.log
#SBATCH -p cpu_medium

### LOAD REQUIRED MODULES
module purge
module load trimgalore/0.5.0 samtools/1.9 bedtools/2.26.0 bowtie2/2.3.4.1 python/cpu/3.7.2 fastqc/0.11.7 igenome/1.0 sratoolkit/2.9.1 

### install more stuff
pip install --user cutadapt
pip install multiqc

### remove files from previous runs
rm *.bam
rm *.sam
rm *.txt
rm *.fq
rm *.zip

rm ./outputs/*

rm ./multiqc_data/*
rmdir multiqc_data

### Download data
fasterq-dump --split-files /gpfs/scratch/bem6982/bioinformatics/assn_02_ngs/fq/fastq_files/SRR7049609.sra
fasterq-dump --split-files /gpfs/scratch/bem6982/bioinformatics/assn_02_ngs/fq/fastq_files/SRR7049610.sra
fasterq-dump --split-files /gpfs/scratch/bem6982/bioinformatics/assn_02_ngs/fq/fastq_files/SRR7049611.sra
fasterq-dump --split-files /gpfs/scratch/bem6982/bioinformatics/assn_02_ngs/fq/fastq_files/SRR7049612.sra
fasterq-dump --split-files /gpfs/scratch/bem6982/bioinformatics/assn_02_ngs/fq/fastq_files/SRR7049615.sra
fasterq-dump --split-files /gpfs/scratch/bem6982/bioinformatics/assn_02_ngs/fq/fastq_files/SRR7049616.sra

### QC TRIMMING OF ADAPTERS AND 3' LOW QUALITY BASES
trim_galore --paired --fastqc --length 30 -q 30 ./fastq_files/SRR7049609.sra_1.fastq ./fastq_files/SRR7049609.sra_2.fastq
trim_galore --paired --fastqc --length 30 -q 30 ./fastq_files/SRR7049611.sra_1.fastq ./fastq_files/SRR7049611.sra_2.fastq
trim_galore --paired --fastqc --length 30 -q 30 ./fastq_files/SRR7049610.sra_1.fastq ./fastq_files/SRR7049610.sra_2.fastq
trim_galore --paired --fastqc --length 30 -q 30 ./fastq_files/SRR7049612.sra_1.fastq ./fastq_files/SRR7049612.sra_2.fastq
trim_galore --paired --fastqc --length 30 -q 30 ./fastq_files/SRR7049615.sra_1.fastq ./fastq_files/SRR7049615.sra_2.fastq
trim_galore --paired --fastqc --length 30 -q 30 ./fastq_files/SRR7049616.sra_1.fastq ./fastq_files/SRR7049616.sra_2.fastq

### FastQC and MultiQC to get html output of trim report
fastqc *.fq ###didn't have to do this last time
multiqc .

### Move htmls to outputs directory
mv *.html ./outputs
