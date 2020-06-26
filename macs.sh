#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 72:00:00
#SBATCH --mem=64g
#SBATCH --mail-type=end
#SBATCH --mail-user=dgoulet@mit.edu

I="KO1_final.bam"
N="KO1"

#load python to cluster to run macs2
module load python/2.7.13

#call peaks in paired-end final bamfile
macs2 callpeak -f BAMPE -t $I -g mm -n $N --keep-dup all --outdir ./macs
