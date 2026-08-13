#!/bin/bash

module load vcftools/0.1.16

vcftools --gzvcf /Data2/amrita/vcf_info/sites_filtered.vcf.gz \
         --thin 10000 \
         --recode \
         --stdout | gzip > /Data2/amrita/vcf_info/thinned_sites.vcf.gz

