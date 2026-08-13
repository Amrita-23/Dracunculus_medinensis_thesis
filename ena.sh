#!/bin/bash
module load ena/
#get fastq files for each accession number for the larvae
for acc_num in $(tail -n +4 GeneticMap_accessions.txt | cut -f 3)
do 
        enaGroupGet -f fastq -d /Data2/amrita/fastq_files "$acc_num"
done
