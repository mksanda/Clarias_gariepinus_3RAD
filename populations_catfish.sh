#!/bin/bash


# ref_map.pl
ref_map.pl -T 8 --samples $HOME/RAD/aligned \
--popmap $HOME/RAD/aligned/popmap/popmap_catfish.tsv \
--out-path $HOME/RAD/aligned/refmap

# Populations analysis
populations --in-path $HOME/RAD/aligned/refmap --out-path $HOME/RAD/aligned/populations/ \
--popmap $HOME/RAD/aligned/popmap/popmap_catfish.tsv \
--threads 8 -r 0.8 -p 9 --min-mac 3 --hwe --write-single-snp --vcf --fstats --genepop
