##Make PCA 


install.packages("SNPRelate", repos="http://R-Forge.R-project.org")
library(SNPRelate)

#Convert vcf to gds for snprelate
snpgdsVCF2GDS(
  "/Users/amritatrehan/documents/final_filtered.vcf.gz",
  "variants.gds",
  method = "biallelic.only"
)
#open gds file
genofile <- snpgdsOpen("variants.gds")

### PLOT THE TERM IBS GRAPH ###
ibs <- snpgdsIBSNum(genofile, autosome.only = FALSE)

ibs.df <- data.frame(
  IBS0 = ibs$ibs0,
  IBS1 = ibs$ibs1,
  IBS2 = ibs$ibs2
)

# Convert to proportions
n <- length(ibs$sample.id)

pairs <- expand.grid(
  sample1 = ibs$sample.id,
  sample2 = ibs$sample.id,
  stringsAsFactors = FALSE
)

pairs$IBS0 <- as.vector(ibs$ibs0)
pairs$IBS1 <- as.vector(ibs$ibs1)
pairs$IBS2 <- as.vector(ibs$ibs2)

total <- pairs$IBS0 + pairs$IBS1 + pairs$IBS2

pairs$IBS0 <- pairs$IBS0 / total
pairs$IBS1 <- pairs$IBS1 / total
pairs$IBS2 <- pairs$IBS2 / total

pairs <- subset(pairs, sample1 != sample2)

library(ggtern)

ggtern(pairs,
       aes(x = IBS0,
           y = IBS1,
           z = IBS2)) +
  geom_point(alpha = 0.25, size = 1,color = 'purple') +
  labs(title = "Pairwise IBS State Ternary Diagram")
  theme_bw()

library(SNPRelate)
snpgdsClose(genofile)




snpset.id <- read.gdsn(index.gdsn(genofile, "snp.id"))
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
##do SNPRELATE calcs


# Estimate IBD coefficients
ibd <- snpgdsIBDMoM(genofile, #sample.id=sample.id, snp.id=snpset.id,
                    maf=0.05, missing.rate=0.05, num.thread=2, autosome.only = FALSE)

# Make a data.frame
ibd.coeff <- snpgdsIBDSelection(ibd)
head(ibd.coeff)

plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1", main="samples (MoM)")
lines(c(0,1), c(1,0), col="red", lty=2)
subset(ibd.coeff, kinship > 0.125)

plot(density(ibd.coeff$kinship))
hist(ibd.coeff$kinship, breaks = 50)

offspring_pairs <- subset(
  ibd.coeff,
  ID1 != "ERS2021560" &
    ID2 != "ERS3380982"
)

hist(offspring_pairs$kinship, breaks = 50)
plot(density(offspring_pairs$kinship))

# Estimating IBD Using Maximum Likelihood Estimation (MLE)
#set.seed(100)
#snp.id <- sample(snpset.id, 1500)  # random 1500 SNPs
ibd <- snpgdsIBDMLE(genofile, #sample.id=sample.id, snp.id=snp.id,
                    maf=0.05, missing.rate=0.05, num.thread=2, autosome.only = FALSE, kinship=TRUE)

# Make a data.frame
ibd.coeff <- snpgdsIBDSelection(ibd)
ibd.coeff$k2=1-(ibd.coeff$k0-ibd.coeff$k1)

install.packages("ggtern")
library(ggtern)
library(ggplot2)
ggtern(data=ibd.coeff,aes(x=k0, y=k1, z=k2))+geom_point()
plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1", main="samples_2")
lines(c(0,1), c(1,0), col="red", lty=2)


packageVersion("ggplot2")
# Incorporate with pedigree information
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

ibd.robust <- snpgdsIBDKING(genofile, sample.id=sample.id,
                            num.thread=2, autosome.only = FALSE)
names(ibd.robust)

# Pairs of individuals
dat <- snpgdsIBDSelection(ibd.robust)

head(dat)

plot(dat$IBS0, dat$kinship,
     xlab = "Proportion of Zero IBS",
     ylab = "Estimated Kinship Coefficient (KING-robust)",
     pch = 1)

plot(dat$IBS0, dat$kinship, xlab="Proportion of Zero IBS",
     ylab="Estimated Kinship Coefficient (KING-robust)")
hist(dat$kinship, breaks = 50)
ibs0 <- ibd.robust$IBS0
kin  <- ibd.robust$kinship
ids  <- ibd.robust$sample.id

dat <- data.frame(
  IBS0 = as.vector(ibs0),
  kinship = as.vector(kin)
)

n <- length(ids)

keep <- upper.tri(ibs0)

dat <- data.frame(
  IBS0 = ibs0[keep],
  kinship = kin[keep]
)

##extract paternal clusters

kin <- ibd.robust$kinship
diag(kin) <- 1

rownames(kin) <- ibd.robust$sample.id
colnames(kin) <- ibd.robust$sample.id

dist_mat <- 1 - kin
hc <- hclust(as.dist(dist_mat), method = "average")

plot(hc, labels = ibd.robust$sample.id,
     main = "Kinship-based clustering")

clusters_2 <- cutree(hc, k = 2)
clusters_3 <- cutree(hc, k = 3)
clusters_4 <- cutree(hc, k = 4)
clusters_5 <- cutree(hc, k = 5)

table(clusters_2)
table(clusters_3)
names(clusters_3[clusters_3 == 3])
names(clusters_4[clusters_4 == 3])
names(clusters_4[clusters_4 == 4])
table(clusters_4)

data.frame(
  sample = names(clusters_3),
  cluster = clusters_3
)

kin <- ibd.robust$kinship

mother <- "ERS2021560"

sort(kin[mother, ], decreasing=TRUE)[1:20]

pedigree <- data.frame(
  individual = ibd.robust$sample.id,
  sire_group = clusters_3
)

write.csv(pedigree, "paternal_families.csv", row.names = FALSE)

#run the pca
pca <- snpgdsPCA(genofile, autosome.only = FALSE)



#variance proportion (%)
head(round(pca$varprop*100, 2))
#13.41  7.30  6.21  4.86  4.65  4.14


#extract the results
sample.id <- pca$sample.id

# make a data.frame
pca.df <- data.frame(
  sample = sample.id,
  PC1 = pca$eigenvect[,1],
  PC2 = pca$eigenvect[,2],
  stringsAsFactors = FALSE)
head(pca.df)

pca.df$cluster <- ifelse(pca.df$PC1 > 0.05, "Group2", "Group1")

table(pca.df$cluster) 

#plot the pca
ggplot(pca.df, aes(PC1, PC2, color= (sample == "ERS2021560"))) +
  geom_point(size = 3) +
  theme_classic()

#plot the pca using father group as colours

pca.df$FatherID <- factor(c("Mother", best_cluster$FatherID))

ggplot(pca.df, aes(PC1, PC2, color = FatherID)) +
  geom_point(size = 3) +
  theme_classic()

#slightly diff plot
plot(pca.df$PC1, pca.df$PC2, xlab="eigenvector 2", ylab="eigenvector 1")


## plot the KINGship analysis but color by PCA 

dat2 <- merge(
  dat,
  pca.df[, c("sample","cluster")],
  by.x = "ID1",
  by.y = "sample",
  all.x = TRUE
)

names(dat2)[names(dat2) == "cluster"] <- "cluster1"

dat2 <- merge(
  dat2,
  pca.df[, c("sample","cluster")],
  by.x = "ID2",
  by.y = "sample",
  all.x = TRUE
)

names(dat2)[names(dat2) == "cluster"] <- "cluster2"

dat2$same_cluster <- dat2$cluster1 == dat2$cluster2


cols <- c(Group1 = "blue",
          Group2 = "red")

plot(dat2$IBS0,
     dat2$kinship,
     col = cols[dat2$cluster1],
     #col = cols[dat2$cluster2],
     xlab = "Proportion of Zero IBS",
     ylab = "Estimated Kinship Coefficient (KING-robust)",
     pch = 19)
head(dat2)

## snpgdsIBSNum will get IBS0,1 and 2


## plot LOD vs number of LGs or number of markers for the biggest 3 LGs
library(ggplot2)
library(tidyr)

# Reproducing your data
data <- tibble(
  LargestLGs = c("Unassigned", "1", "2", "3"),
  Lod_5 = c(34, 7377, 379, 0),
  Lod_6 = c(58, 7251, 379, 102),
  Lod_7 = c(93, 7195, 378, 102),
  Lod_8 = c(153, 7043, 377, 101),
  Lod_9 = c(250, 6648, 336, 101),
  Lod_10 = c(362, 4780, 1438, 296),
  Lod_11 = c(496, 1531, 1159, 954),
  Lod_12 = c(678, 1119, 995, 878)
)

to_plot <- pivot_longer(data, cols = starts_with("Lod"), names_to = "Zone", values_to = "Count_of_Markers")

ggplot(to_plot, aes(x = LargestLGs, y = Count_of_Markers, fill = Zone)) +
  geom_text(
    aes(label = Count_of_Markers),
    position = position_dodge(0.9),
    size = 2.5,
    hjust = -0.2,
    angle = 90
  ) +
  geom_bar(position="dodge", stat="identity") + labs(title="LOD score affect on Linkage Group Size", y="Number of Markers", x="Linkage Groups")

## plot datafilter vs number of LGs or number of markers for the biggest 3 LGs
library(ggplot2)
library(tidyr)

# Reproducing your data
data <- tibble(
  LargestLGs = c("Unassigned", "1", "2", "3"),
  "0.1" = c(3786, 900, 900, 626),
  "0.01" = c(2578, 981, 928, 895),
  "0.001" = c(1755, 984, 931, 915),
  "0.0001" = c(872, 1015, 981, 935),
  "0.00001" = c(405, 1026, 983, 983),
  "0.000001" = c(405, 947, 919, 828),
  "0.0000001" = c(38, 1395, 1415, 1361)
)

to_plot <- pivot_longer(data, cols = starts_with("0."), names_to = "Zone", values_to = "Count_of_Markers")

ggplot(to_plot, aes(x = LargestLGs, y = Count_of_Markers, fill = Zone)) +
  geom_text(
    aes(label = Count_of_Markers),
    position = position_dodge(0.9),
    size = 2.5,
    hjust = -0.2,
    angle = 90
  ) +
  geom_bar(position="dodge", stat="identity") + labs(title="Data Tolerance affect on Linkage Group Size", y="Number of Markers", x="Linkage Groups")
