#!/bin/bash
set -euo pipefail

module load slurmzy/
module load samtools/1.17
module load gatk/4.1.3.0

#get reference
REF=/Data2/amrita/dracunculus_medinensis_v3.fa

# combine all ERS gVCFs into one cohort gVCF
variants=""

for gvcf in /Data2/amrita/vcf/*.g.vcf.gz; do
    [[ -f "$gvcf" ]] || continue
    variants="$variants --variant $gvcf"
done

gatk CombineGVCFs \
    -R "$REF" \
    $variants \
    -O /Data2/amrita/vcf/cohort.g.vcf.gz


# convert final cohort gVCF into VCF
gatk GenotypeGVCFs \
    -R "$REF" \
    -V /Data2/amrita/vcf/cohort.g.vcf.gz \
    -O /Data2/amrita/vcf/cohort.vcf.gz
