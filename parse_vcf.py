#open vcf 
from cyvcf2 import VCF, Writer
from collections import defaultdict

vcf = VCF("/Users/amritatrehan/Documents/final_filtered.vcf.gz")
w = Writer("sites_filtered.vcf.gz", vcf)
sample_name = vcf.samples[0] #this is the mom! 
#print(sample_name)
#make dict of mom het sites
mom_hets={}
contig_sites = defaultdict(list)
#get heterozygeous sites for this sample

three_to_one_sites = {}
mother = 0 

for variant in vcf:

    if variant.gt_types[mother] != 1:
        continue

    key = f"{variant.CHROM}:{variant.POS}"

    #quality control the mother so qual is greater than 30 
    gq = variant.format("GQ")
    if gq is not None:
        m_gq = gq[mother][0]
        if m_gq < 25:
            continue

    mom_hets[key] = {
        "ref": variant.REF,
        "alt": ",".join(variant.ALT)
    }
#print(len(mom_hets))
    alt = 0
    total = 0
    # all of the samples except the mom
    for i in range(1, len(vcf.samples)):
        #genotypes of each sample
        gt = variant.genotypes[i]
        # the two alleles
        a1, a2 = gt[:2]

        #make sure all of the alleles are real
        if a1 < 0 or a2 < 0:
            continue
        
        #total number of ALT alleles across all offspring
        alt += (a1 == 1) + (a2 == 1)
        #count total alleles
        total += 2
    #if the genotype is not valid skip it 
    if total == 0:
        continue

    # calculate al allele frequency 
    freq = alt / total

    # keep ~3:1 sites

    if 0.60 <= freq <= 0.85 or 0.15 <= freq <= 0.40:
        w.write_record(variant)
        contig_sites[variant.CHROM].append({
            "pos": variant.POS,
            "freq": freq
        })
        
        #count number of sites 
        three_to_one_sites[f"{variant.CHROM}:{variant.POS}"] = {
            "ref": variant.REF,
            "alt": ",".join(variant.ALT),
            "freq": freq
        }

print("3:1 candidate sites:", len(three_to_one_sites)) #99805 sites 
for chrom, positions in contig_sites.items():
    print(chrom, len(positions))

# get list of sites for vcf positions to use for lepmap3
with open("/Users/amritatrehan/Documents/sites.txt", "w") as out:
    for site in three_to_one_sites:
        chrom, pos = site.split(":")
        out.write(f"{chrom}\t{pos}\n")


w.close()
vcf.close()
