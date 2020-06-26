#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 24:00:00
#SBATCH --mem=64g
#SBATCH --mail-type=end
#SBATCH --mail-user=dgoulet@mit.edu

#load fastqc package to cluster
module add fastq

#run fastqc
fastqc /home/dgoulet/data/ATAC/PHF6/test1.fastq.gz --outdir=./reports
