#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 72:00:00
#SBATCH --mem=64g
#SBATCH --mail-type=end
#SBATCH --mail-user=dgoulet@mit.edu

#Designate filepath for input and output files
I="KO1_filt_dedup_mito.bam"
O="KO1_filt_dedup_mito_tn5.bam"
F="KO1_final.bam"

#load python deeptools and samtools packages to cluster
module load python/2.7.13
module load deeptools
module load samtools

#Shift reads to correct for Tn5 integration
alignmentSieve -p 8 --ATACshift -b $I -o $O

#Sort corrected bamfile
samtools sort -o $F -@ 8 $O

#Index sorted and corrected bamfile
samtools index -@ 8 $F
