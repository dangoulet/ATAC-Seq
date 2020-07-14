#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 72:00:00
#SBATCH --mem=64g
#SBATCH --mail-type=end
#SBATCH --mail-user=kerberos@mit.edu

#load packages to cluster
module load fastqc
module load bedtools
module load python/2.7.13
module load deeptools
module load samtools
module load bowtie2/2.3.5.1
module load picard
module load r/3.6.0

#designate variables for files
R1="./raw/KO1_mix_S4_R1_001.trim.fastq.gz"
R2="./raw/KO1_mix_S4_R2_001.trim.fastq.gz"
SAMPLE="S1"
idx="/home/dgoulet/data/Genomes/mm10/Sequence/Bowtie2Index/genome"
M="/home/dgoulet/data/Blacklist/mm10_chrM.bed"
bl="home/dgoulet/data/Blacklist/mm10blacklist.bed"

#make reports directory if it does not already exist
mkdir -p reports

#run fastqc
fastqc $R1 $R2 -o ./reports/

#Use very sensitive alignment then sort reads and output to bamfile
bowtie2 --very-sensitive -p 8 -x $idx -1 $R1 -2 $R2 | samtools view  -u  - | samtools sort  - > ${SAMPLE}.bam

#filter for mapped reads with unambiguous alignment
samtools view  -b  -@ 8  -q 10  -F 4  ${SAMPLE}.bam  >  ${SAMPLE}_filt.bam

#Output bamfile stats to textfile to check filtering
samtools flagstat ${SAMPLE}_filt.bam > ./reports/${SAMPLE}_filter_report.txt

#filter out duplicate reads
java -jar $PICARD MarkDuplicates I=${SAMPLE}_filt.bam O=${SAMPLE}_filt_dedup.bam M=./reports/${SAMPLE}_dups.txt REMOVE_DUPLICATES=true

#export bamfile stats to textfile
samtools flagstat ${SAMPLE}_filt_dedup.bam > ./reports/${SAMPLE}_dedup_stats.txt

#filter chrM and chrUn from bamfile
samtools view  -@ 8  -h  ${SAMPLE}_filt_dedup.bam | awk '($3 != "chrM" && $3 != "chrUn")' | samtools view  -@ 8  -hb  - > ${SAMPLE}_filt_dedup_mito.bam

#index bamfile
samtools index  -@ 8  -b  ${SAMPLE}_filt_dedup_mito.bam

#export chrM coverage to textfile
samtools bedcov $M ${SAMPLE}_filt_dedup_mito.bam > ./reports/${SAMPLE}_filt_dedup_mito_stats.txt

#export ATAC insert size to textfile and histogram
java -jar $PICARD CollectInsertSizeMetrics I=${SAMPLE}_filt_dedup_mito.bam H=./reports/${SAMPLE}_insert_hist.pdf O=./reports/${SAMPLE}_insert_metrics.txt

#Shift reads to correct for Tn5 integration
alignmentSieve -p 8 --ATACshift -b ${SAMPLE}_filt_dedup_mito.bam -o ${SAMPLE}_filt_dedup_mito_tn5.bam

#Sort corrected bamfile
samtools sort  -o ${SAMPLE}_final.bam  -@ 8  ${SAMPLE}_filt_dedup_mito_tn5.bam

#Index sorted and corrected bamfile
samtools index  -@ 8  ${SAMPLE}_final.bam

#Modify the timestamp to prevent downstream errors
touch -m ${SAMPLE}_final.bam.bai

#call peaks in paired-end final bamfile
macs2 callpeak -f BAMPE -t ${SAMPLE}_final.bam -g mm -n ${SAMPLE} --keep-dup all --outdir ./macs

bedtools intersect -nonamecheck -v -a ./macs/${SAMPLE}_peaks.narrowPeak -b $bl > ./macs/${SAMPLE}_peaks_bl.narrowPeak
