#!/bin/bash

mkdir -p /Data2/amrita/bwa_files

REF=/Data2/amrita/dracunculus_medinensis_v3.fa
___
R1=${1}
        #make read 2 matching file
        R2=${R1/_1.fastq.gz/_2.fastq.gz}

        #get sample name
        SAMPLE=$(basename "$R1" _1.fastq.gz)

        echo "Processing $SAMPLE"

        #do the actual alignment and convert sam files to bam files then sort them
        bwa mem -M -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
         "$REF" "$R1" "$R2" | samtools sort -o /Data2/amrita/bwa_files/${SAMPLE}.sorted.bam

