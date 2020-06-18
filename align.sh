#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 72:00:00
#SBATCH --mem=64g
#SBATCH --mail-type=end
#SBATCH --mail-user=dgoulet@mit.edu

#Designate filepath for input files and index
R1="test1.fastq"
R2="test2.fastq"
idx="/home/dgoulet/data/Genomes/mm10/Sequence/Bowtie2Index/genome"

#load bowtie2 and samtools packages to cluster
module add bowtie2/2.3.5.1

module add samtools/1.5

#Use very sensitive alignment then sort reads and output to bamfile
bowtie2 --very-sensitive -p 8 -x $idx -1 $R1 -2 $R2 | samtools view -u - | samtools sort - > test.bam
