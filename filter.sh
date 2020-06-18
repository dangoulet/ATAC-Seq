#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 72:00:00
#SBATCH --mem=64g
#SBATCH --mail-type=end
#SBATCH --mail-user=dgoulet@mit.edu

#Designate filepath for input and output files
I="KO1.bam"
O="KO1_filt.bam"

#load samtools package to cluster
module load samtools/1.5

#filter for mapped reads with unambiguous alignment
samtools view -b  -q 10  -F 4  $I  >  $O

#Output bamfile stats to textfile to check filtering
samtools flagstat $O > KO1_filter_report.txt
