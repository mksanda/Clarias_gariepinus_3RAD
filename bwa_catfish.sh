#!/bin/bash

set -x
#index the genome
bwa index $HOME/RAD/ref_db/catfish/clarias.fa 

#change directory to access fastq files 
cd $HOME/RAD/fastq/

# get list of all file names, generated separately
SAMPLELIST=`ls -1 *.1.fa.gz | sed -E 's/.1.fa.gz//'`

# align each pair of sample reads to the reference genome
for SAMPLE in `cat $SAMPLELIST`; do

bwa mem -t 8 $HOME/RAD/ref_db/catfish/clarias \
/data/ref_db/catfish/fq_files/${SAMPLE}.1.fq.gz \
/data/ref_db/catfish/fq_files/${SAMPLE}.2.fq.gz \
| samtools view -h -b | \
 samtools sort --thread 8 > $HOME/RAD/aligned/${SAMPLE}.bam

done

