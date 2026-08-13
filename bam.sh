#!/bin/bash
module load picard/
module load gatk/
#make new directory to store bam/sam files

mkdir -p /Data2/amrita/sample_bams
REF=/Data2/amrita/dracunculus_medinensis_v3.fa
Bam_dir=/Data2/amrita/bwa_files/
new_bamdir=/Data2/amrita/sample_bams

#get the textfile
find /home/amrita/sample_names/ERS2021560.txt | while read -r textfile
#find /home/amrita/sample_names -name "*.txt" | while read -r textfile
do
        SAMPLE=$(basename "$textfile" .txt)
        slurmzy run 32 "$SAMPLE" "bash new_bam.sh '$textfile'"
done
