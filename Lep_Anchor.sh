#!/bin/bash
OUTDIR=/Data2/amrita/mapping_files
BIN=/Data2/amrita/vcf_info/Lep_anchor/bin
AWK=/Data2/amrita/vcf_info/Lep_anchor/awk

paste snps.txt map11_js.txt | awk 'NR>1' > cleanMap.input

java -cp /Data2/amrita/vcf_info/Lep_anchor/bin \
CleanMap \
map=contig_pos_map_sorted.txt \
> ${OUTDIR}/map.clean

java -cp /Data2/amrita/vcf_info/Lep_anchor/bin \
Map2Bed \
map= ${OUTDIR}/map.clean \
contigLength=/Data2/amrita/dracunculus_medinensis_v3.sizes \
>  ${OUTDIR}/map.bed

awk '$5==1'  ${OUTDIR}/map.bed >  ${OUTDIR}/chr1.bed
awk '$5==2'  ${OUTDIR}/map.bed >  ${OUTDIR}/chr2.bed
awk '$5==3'  ${OUTDIR}/map.bed >  ${OUTDIR}/chr3.bed


for CHROM in {1..3}
do
        java -cp /Data2/amrita/vcf_info/Lep_anchor/bin \
        PlaceAndOrientContigs \
        map= ${OUTDIR}/LG${CHROM}_place.map \
        bed= ${OUTDIR}/chr${CHROM}.bed \
        >  ${OUTDIR}/LG${CHROM}.la \
        2>  ${OUTDIR}/LG${CHROM}.la.err
done

for CHROM in {1..3}
do
        awk -vlg=${CHROM} \
        -f /Data2/amrita/vcf_info/Lep_anchor/awk/makeagp_full2.awk \
         ${OUTDIR}/LG${CHROM}.la \
        >  ${OUTDIR}/LG${CHROM}.agp
done
for CHROM in {1..3}
do
        awk -vn=X \
        '(NR==FNR){map[NR-1]=$0}
        (NR!=FNR){$1=map[$1] "\t" n; print}' \
         ${OUTDIR}/snps.txt  ${OUTDIR}/intervals_LG${CHROM}.int \
        >  ${OUTDIR}/order${CHROM}_i.input
done

for CHROM in {1..3}
do
        awk \
        -f ${AWK}/liftover.awk \
         ${OUTDIR}/LG${CHROM}.agp \
         ${OUTDIR}/order${CHROM}_i.input | \
        sort -V | \
        grep CHR \
        >  ${OUTDIR}/order${CHROM}.liftover
done

for CHROM in {1..3}
do

        awk -vinverse=1 \
        -f ${AWK}/liftover.awk \
         ${OUTDIR}/LG${CHROM}.agp \
         ${OUTDIR}/order${CHROM}.liftover | \
        awk '(NR==FNR){
                m[$1"\t"($2+0)]=NR-1
        }
                (NR!=FNR){
                         print m[$1"\t"($2+0)]
        }' \
        snps.txt - \
        >  ${OUTDIR}/order${CHROM}.phys

done
