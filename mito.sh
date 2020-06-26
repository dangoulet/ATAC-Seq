#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 72:00:00
#SBATCH --mem=64g
#SBATCH --mail-type=end
#SBATCH --mail-user=dgoulet@mit.edu

#designate filepath of input, output, as well as mm10 chrM bedfile
I="KO1_filt_dedup.bam"
O="KO1_filt_dedup_mito.bam"
M="/home/dgoulet/data/Blacklist/mm10_chrM.bed"

#load picard and samtools packages to cluster
module load samtools

#filter chrM and chrUn from bamfile
samtools view -h $I | awk '($3 != "chrM" && $3 != "chrUn")' | samtools view -hb - > $O

#index bamfile
samtools index -b $O

#export chrM coverage to textfile
samtools bedcov $M $O > ./reports/KO1_filt_dedup_mito_stats.txt
