#!/bin/bash

# Set values for input parameters
raw_data=$HOME/RAD/raw/
output_folder=$HOME/RAD/fastq/
library_info=$HOME/RAD/barcodes_3rad.txt

##helpful link https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php
process_radtags -P -c -q -r -p ${raw_data} -o ${output_folder} -b ${library_info} --inline_inline --renz_1 mspI --renz_2 bamHI -i gzfastq


