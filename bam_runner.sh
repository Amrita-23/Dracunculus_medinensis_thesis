#!/bin/bash

module load slurmzy/
module load gatk/
module load samtools/

fastq=/Data2/amrita/fastq_files
for f in "$fastq"/*; do
        sample=$(basename "$f")
        {
                touch "$sample".txt

                printf "%s\t" "$sample"

                for sub in "$f"/*/; do
                        printf "%s " "$(basename "$sub")"
                done
                printf "\n" 
        } > "${sample}.txt"


done
