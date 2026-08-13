#!/bin/bash
module load gatk/4.1.3.0

gatk SelectVariants \
  -R /Data2/amrita/dracunculus_medinensis_v3.fa \
  -V /Data2/amrita/vcf_info/snps_filtered.vcf.gz \
  --exclude-filtered true \
  -O /Data2/amrita/vcf_info/pass_snps_filtered.vcf.gz
