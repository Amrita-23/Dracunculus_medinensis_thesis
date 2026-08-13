#LepMap3 code used to create female recombination map run in Linux:

#(For male recombination map informativeMask was set to 13, and the pedigree built for male #recombination was used)

LEPMAP=/Data2/amrita/vcf_info/LepMap3/bin
DATA=/Data2/amrita/mapping_files/data_f_min5.call
OUTDIR=/Data2/amrita/mapping_files

java -cp ${LEPMAP} ParentCall2 ignoreParentOrder=1 \ 
data=/Data2/amrita/vcf_info/pedigree_fixed.jac.txt \ 
vcfFile=/Data2/amrita/vcf_info/thinned_sites.vcf >/Data2/amrita/mapping_files/data.call

java -cp /Data2/amrita/vcf_info/LepMap3/bin Filtering2 
data=${OUTDIR}/data.call
removeNonInformative=1 \ 
dataTolerance=0.0000001 \
 >${OUTDIR}/data_f.call \

# get rid of scaffolds with less than five markers 

awk 'NR>7 {print $1}' ${OUTDIR}/data_f.call | sort | uniq -c | awk '$1>=5 {print $2}' > ${OUTDIR}/keep_scaffolds.txt
awk ' FNR==NR {keep[$1]=1; next} FNR<=7 {print; next} $1 in keep '${OUTDIR}/ keep_scaffolds.txt ${OUTDIR}/data_f.call > ${OUTDIR}/data_f_min5.call
# lines below run in separate scripts # 24-49
    # Separate chromosomes
java -cp ${LEPMAP} SeparateChromosomes2 \
        data=${OUTDIR}/data_f_min5.call \
        lodLimit=11 \
        distortionLod=1 \
        > ${OUTDIR}/map11_0.0000001.txt

    # Join singles
java -cp ${LEPMAP} JoinSingles2All \
        map=${OUTDIR}/map11_0.0000001.txt \
        data=${OUTDIR}/data_f_min5.call \
        lodLimit=2 \
        > ${OUTDIR}/map11_js_0.0000001.txt
for CHROM in 1 2 3
do
    # Order LG1
        java -cp ${LEPMAP} OrderMarkers2 \
                data=${OUTDIR}/data_f_min5.call \
                map=${OUTDIR}/map11_js_0.0000001.txt \
                chromosome=${CHROM} \
                usePhysical=1 \
                informativeMask=2 \
                hyperPhaser=1 \
                proximityScale=50 \
                > ${OUTDIR}/order${CHROM}_usephys_0.0000001.txt
done
##run in its own script after the first two were run

cat data_f_min5.call|cut -f 1,2|awk '(NR>=7)' >snps_new.txt
for CHROM in 1 2 3
do
    # Order LG1
        java -cp ${LEPMAP} OrderMarkers2 \
                data=${DATA} \
                map=${OUTDIR}/map11_js_0.0000001.txt \
                chromosome=${CHROM} \
                informativeMask=2 \
                hyperPhaser=1 \
                improveOrder=0 \
                proximityScale=50 \
                computeLODScores=${OUTDIR}/0.0000001_LOD_scores${CHROM}.txt \
                evaluateOrder=${OUTDIR}/order${CHROM}_usephys_0.0000001.txt \
                calculateIntervals=${OUTDIR}/0.0000001_intervals_LG${CHROM}.int \
                >  ${OUTDIR}/order${CHROM}_after.txt
        #create map
        awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' ${OUTDIR}/snps_new.txt \
        ${OUTDIR}/order${CHROM}_after.txt \
        >${OUTDIR}/order${CHROM}_eval_0.0000001.mapped
done
