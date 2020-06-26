#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 72:00:00
#SBATCH --mem=64g
#SBATCH --mail-type=end
#SBATCH --mail-user=dgoulet@mit.edu

#Designate filepath for input and output files
I="test.bam"
O="test_dedup.bam"

#load picard and samtools packages to cluster
module load picard
module load samtools/1.5

#filter out duplicate reads
java -jar $PICARD MarkDuplicates I=$I O=$O M=dups.txt REMOVE_DUPLICATES=true

#export bamfile stats to textfile
samtools flagstat $O > ./reports/test_dedup_stats.txt
