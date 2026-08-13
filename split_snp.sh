#!/bin/bash
module load gatk/4.1.3.0

gatk SelectVariants \
    -R /Data2/amrita/dracunculus_medinensis_v3.fa \
    -V /Data2/amrita/vcf/cohort.vcf.gz \
    --select-type-to-include SNP \
    -O /Data2/amrita/vcf/snp_cohort.vcf.gz

