#!/bin/bash
module load gatk/4.1.3.0

gatk VariantFiltration \
        -R /Data2/amrita/dracunculus_medinensis_v3.fa \
        -V /Data2/amrita/vcf/snp_cohort.vcf.gz \
        -filter "QD<9.0" --filter-name LowQD \
        -filter "QUAL<30.0" --filter-name LowQUAL \
        -filter "SOR>2.83" --filter-name HighSOR \
        -filter "FS>18.5" --filter-name HighFS \
        -filter "MQ<48.5" --filter-name LowMQ \
        -filter "MQRankSum<-1.15" --filter-name LowMQRankSum \
        -filter "ReadPosRankSum<-2.55" --filter-name LowReadPosRankSum \
        -O /Data2/amrita/vcf_info/snps_filtered.vcf.gz
