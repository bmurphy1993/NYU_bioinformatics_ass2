#!/bin/bash
#SBATCH --job-name=btw_sam_11
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bem6982@nyu.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/scratch/bem6982/bioinformatics/assn_02_ngs/logs/NGS1_%j.log
#SBATCH -p cpu_medium


### Modules
module purge
module load trimgalore/0.5.0 samtools/1.9 bedtools/2.26.0 bowtie2/2.3.4.1 python/cpu/3.7.2 fastqc/0.11.7 igenome/1.0

### bowtie2
bowtie2 -p 16 -x /gpfs/share/apps/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome -1 SRR7049611.sra_1_val_1.fq -2 SRR7049611.sra_2_val_2.fq -S SRR7049611.mapping.sam

### SAMTOOLS PROCESSING
samtools view -b -o SRR7049611.mapping.bam SRR7049611.mapping.sam
samtools sort -o SRR7049611.mapping.sorted.bam SRR7049611.mapping.bam
samtools index SRR7049611.mapping.sorted.bam
cp SRR7049611.mapping.sorted.bam ./outputs
cp SRR7049611.mapping.sorted.bam.bai ./outputs

### STRAND SEPERATION
samtools view -b -f 83 SRR7049611.mapping.sorted.bam > fwd1_11.bam
samtools view -b -f 163 SRR7049611.mapping.sorted.bam > fwd2_11.bam
samtools index fwd1_11.bam
samtools index fwd2_11.bam
samtools merge -f fwd_11.bam fwd1_11.bam fwd2_11.bam
samtools index fwd_11.bam
cp fwd_11.bam ./outputs
cp fwd_11.bam.bai ./outputs

samtools view -b -f 99 SRR7049611.mapping.sorted.bam > rev1_11.bam
samtools view -b -f 147 SRR7049611.mapping.sorted.bam > rev2_11.bam
samtools index rev1_11.bam
samtools index rev2_11.bam
samtools merge -f rev_11.bam rev1_11.bam rev2_11.bam
samtools index rev_11.bam
cp rev_11.bam ./outputs
cp rev_11.bam.bai ./outputs

### BEDGRAPH FILES
samtools view -b fwd_11.bam | genomeCoverageBed -ibam stdin -bg -split > ./outputs/SRR7049611-fwd.bedgraph
samtools view -b rev_11.bam | genomeCoverageBed -ibam stdin -bg -split > ./outputs/SRR7049611-rev.bedgraph
samtools view -b SRR7049611.mapping.sorted.bam | genomeCoverageBed -ibam stdin -bg -split > ./outputs/SRR7049611.bedgraph

### BED12 FILES
bamToBed -bed12 -i fwd_11.bam > ./outputs/SRR7049611-fwd.bed
bamToBed -bed12 -i rev_11.bam > ./outputs/SRR7049611-rev.bed
bamToBed -bed12 -i SRR7049611.mapping.sorted.bam > ./outputs/SRR7049611.bed
