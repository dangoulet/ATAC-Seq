# ATAC-Seq
Pipeline for Local ATAC-Seq Processing on Luria Compute Cluster
Dan Goulet Hemann Lab MIT 2020

This ATAC-Seq Pipeline is a series of shell scripts written to process ATAC-Seq data in a stepwise manner:

1. fastq.sh - This script will quantify the quality of reads in your fastq files.

2. align.sh - This script aligns your reads to the genome index located at the filepath designated by the idx variable.

3. filter.sh - This script filters for mapped reads with unambiguous alignments, and records the output file stats to a textfile.

4. dedup.sh - This script filters duplicate reads from the aligned bamfile and records the output file stats to a textfile.

5. mito.sh - This script filters any reads aligned to the mitochondrial chromosome or an unknown chromosome, then indexes the output and records the number of reads mapping to the mitochondrial chromosome in the output file.

6. Tn5.sh - This script corrects for Tn5 integration, shifting reads +5bp on the positive strand, and -4bp on the negative strand, then indexes the output file.

7. macs.sh - This script calls peaks in the final aligned, filtered, sorted, indexed, and shifted bamfile.

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
