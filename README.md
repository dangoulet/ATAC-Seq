# ATAC-Seq
Pipeline for Local ATAC-Seq Processing on Luria Compute Cluster
Dan Goulet Hemann Lab MIT 2020

This ATAC-Seq Pipeline is a series of scripts written to process ATAC-Seq data in a stepwise manner:

1. Designate variables - You must designate the filepath for your forward read (R1) and reverse read (R2) (e.g. "./raw/KO1_mix_S4_R1_001.trim.fastq.gz"). Then, designate a unique sample prefix (e.g. "S1") that will distinguish this sample. This prefix will be added to filenames at each step, so multiple samples from the sample folder can be processed in parallel without error. Finally, designate the filepath for the Bowtie2 index (e.g. "/home/dgoulet/data/Genomes/mm10/Sequence/Bowtie2Index/genome") used for sequence alignment, and the filepath for the bedfile containing the mitochondrial chromosome location (e.g. "/home/dgoulet/data/Blacklist/mm10_chrM.bed")

2. Check read quality - Fastqc will quantify the quality of reads in your fastq files and record the output in an html file in the reports folder

3. Align reads - Bowtie2 aligns your reads to the genome index located at the filepath designated by the idx variable. This script is written for multi-threaded use with 8 processors. If you change the number of processors requested in your SLURM script (e.g. #SBATCH -n 8), make sure you change the -p option to reflect the new number of processors. --very-sensitive option aligns only the highest quality reads. If the quality of the sequencing run is low, you can use an alternative, found in the bowtie2 manual.

3. Filter reads - This step uses samtools view to identify mapped reads with unambiguous alignments (-F 4 = mapped reads, -q 10 = unabiguous alignment), and records the output file stats to a textfile.

4. Remove duplicates - MarkDuplicates filters duplicate reads from the aligned bamfile and records the output file stats to a textfile.

5. Filter mitochondrial reads - This step uses samtools view to identify mapped reads, then runs a simple bash command to filter any reads aligned to the mitochondrial chromosome or an unknown chromosome. The output is then indexed and the number of reads mapping to the mitochondrial chromosome is recorded in an output file.

6. Record Tn5 insert size - This step uses CollectInsertSizeMetrics to quatify the size of the Tn5 insert and record the output in a text file and histogram.

6. Tn5 Shift - This step uses deeptools alignmentSieve package to correct for Tn5 integration, shifting reads +5bp on the positive strand, and -4bp on the negative strand. The output is then indexed using the samtools index command.

7. Modify index timestamp - the touch -m command modifies the timestamp of the bamfile index, which can throw an error in downstream processing if the timestamp is not updated.

8. Call peaks - This step uses MACS2 to calls peaks in the final aligned, filtered, sorted, indexed, and shifted bamfile.

Package options used in each script are explained in the resources below.

Resources

bowtie2
http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

samtools view
http://www.htslib.org/doc/samtools-view.html

samtools flagstat
http://www.htslib.org/doc/samtools-flagstat.html

samtools coverage
http://www.htslib.org/doc/samtools-coverage.html

samtools sort
http://www.htslib.org/doc/samtools-sort.html

bedtools intersect
https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html

picard CollectInsertSizeMetrics
https://broadinstitute.github.io/picard/command-line-overview.html#CollectInsertSizeMetrics

picard MarkDuplicates
https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates

macs2 callpeak
https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html

Harvard Informatics Pipeline
https://informatics.fas.harvard.edu/atac-seq-guidelines-old-version.html
