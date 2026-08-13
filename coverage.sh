#!/bin/bash#

module load slurmzy/0.0.1
module load bedtools/2.31

mkdir -p /Data2/amrita/bed_files

bedtools makewindows \
  -g /Data2/amrita/dracunculus_medinensis_v3.fa.fai \
  -w 10000 \
  > /Data2/amrita/bed_files/genome_10kb.bed

find /Data2/amrita/sample_bams -name "*.bam" | while read -r R1
do
        SAMPLE=$(basename "$R1" .bam)

slurmzy run 32 "$SAMPLE" "bash /home/amrita/coverage_runner.sh \"$R1\""
done
