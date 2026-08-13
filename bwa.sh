#!/bin/bash

#align fastq reads to the reference genome

#get reference
REF=/Data2/amrita/dracunculus_medinensis_v3.fa
#index ref
bwa index "$REF"

#make new directory to store bam/sam files
mkdir -p /Data2/amrita/bwa_files

#load all of the stuff
load module bwa/0.7.17
load module slurmzy
load module samtools/1.17

#get the fastq files
#find /Data2/amrita/fastq_files -name "*_1.fastq.gz" | while read -r R1
find /Data2/amrita/fastq_files -name "ERR3083828_1.fastq.gz" | while read -r R1
do
        SAMPLE=$(basename "$R1" _1.fastq.gz)
        slurmzy run 32 "$SAMPLE" "bash bwa_runner.sh $R1"
done
