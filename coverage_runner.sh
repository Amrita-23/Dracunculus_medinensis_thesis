#!/bin/bash

module load bedtools/2.31

#for bam in "$@"; do
for bam in /Data2/amrita/sample_bams/*.bam; do
        sample=$(basename "$bam" .bam)
        bedtools coverage \
                -a /Data2/amrita/bed_files/genome_10kb.bed \
                -b "$bam" \
                -sorted \
                -mean \
                > /Data2/amrita/bed_files/"${sample}.10kb_cov.txt"
done
