#!/bin/bash 

module load picard/2.27.5
module load gatk/4.1.3.0
#make new directory to store bam/sam files

mkdir -p /Data2/amrita/reheadered_bams/
mkdir -p /Data2/amrita/sample_bams/

REF=/Data2/amrita/dracunculus_medinensis_v3.fa
Bam_dir=/Data2/amrita/bwa_files/
re_bams=/Data2/amrita/reheadered_bams/
new_bamdir=/Data2/amrita/sample_bams/
textfile=$1



while IFS=$'\t' read -r column1 column2; do
#you can then set IFS to ' ' and loop through column2 entries:
#echo "c1:"$column1
#echo "c2:"$column2
        merge_args=()
        BAMS=()
        read -r -a BAMS <<< "$column2"

        for i in ${BAMS[@]}; do
                echo "i:"$i
                #do stuff per BAM
                picard MarkDuplicates \
                        -I ${Bam_dir}${i}.sorted.bam \
                        -O ${Bam_dir}${i}_marked_duplicates.bam \
                        -M marked.txt
                picard AddOrReplaceReadGroups \
                        -I ${Bam_dir}${i}_marked_duplicates.bam \
                        -O ${re_bams}/${column1}_${i}.bam \
                        -RGID $i \
                        -RGLB $i \
                        -RGPL ILLUMINA \
                        -RGPU 1 \
                        -RGSM $column1
                merge_args+=(-I "${re_bams}/${column1}_${i}.bam")

        done
        picard MergeSamFiles \
                "${merge_args[@]}" \
                -O "${new_bamdir}/${column1}.bam"
        samtools index "${new_bamdir}/${column1}.bam"


# do stuff per sample
        gatk HaplotypeCaller \
                -R "$REF" \
                -I "${new_bamdir}${column1}.bam" \
                -O "/Data2/amrita/vcf/${column1}.g.vcf.gz" \
                -ERC GVCF

        gatk IndexFeatureFile \
                -F "/Data2/amrita/vcf/${column1}.g.vcf.gz"



done < "$textfile"  #put file name here
