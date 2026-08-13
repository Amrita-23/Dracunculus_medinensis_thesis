#!/bin/bash
module load tabix/1.3.1
module load vcftools/0.1.16
#filter out missing to samples with more than 5% missing are filtered out and the missing indiv are also filtered from the snps that already passed the hard filtering

vcftools \
  --gzvcf /Data2/amrita/vcf_info/pass_snps_filtered.vcf.gz \
  --remove missing_samples.txt \
  --max-missing 0.95 \
  --recode --stdout \
  |bgzip -c > /Data2/amrita/vcf_info/final_filtered.vcf.gz

tabix -p vcf /Data2/amrita/vcf_info/final_filtered.vcf.gz

