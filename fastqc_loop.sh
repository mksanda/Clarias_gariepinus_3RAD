#!/usr/bin/bash


# Set values for input parameters
raw=/mnt/d/catfish_process_radtags/
fastqc_F=/mnt/d/catfish_process_radtags/fastqc_ff/
fastqc_R=/mnt/d/catfish_process_radtags/fastqc_rr/
multiqc_F=/mnt/d/catfish_process_radtags/multiQC_ff/
multiqc_R=/mnt/d/catfish_process_radtags/multiQC_rr/

#forward read
for f in ${raw}*1.fq.gz;
do fastqc $f -o ${fastqc_F}
done


#reverse read
for f in ${raw}*2.fq.gz;
do fastqc $f -o ${fastqc_R}
done



