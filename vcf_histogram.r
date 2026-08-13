annotation_table <- read.table("/Users/amritatrehan/documents/annotation_table.txt", col.names=c("CHROM", "POS", "REF", "ALT", "QUAL", "MQ", "MQRankSum", "ReadPosRankSum", "QD", "SOR", "FS"), sep="\t")
missing_indiv <- read.table("/Users/amritatrehan/documents/missing_indv.imiss",  header=TRUE)
library(ggplot2)
missing_site <- read.table("/Users/amritatrehan/documents/missing_site.lmiss", header=TRUE)


#make everything into a numeric
annotation_table$QUAL <- as.numeric(annotation_table$QUAL)
annotation_table$QD <- as.numeric(annotation_table$QD)
annotation_table$FS <- as.numeric(annotation_table$FS)
annotation_table$SOR <- as.numeric(annotation_table$SOR)
annotation_table$MQ <- as.numeric(annotation_table$MQ)
annotation_table$MQRankSum <- as.numeric(annotation_table$MQRankSum)
annotation_table$ReadPosRankSum <- as.numeric(annotation_table$ReadPosRankSum)

#plot histograms for each variable dist
##QUAL 

ggplot(annotation_table, aes(x = (QUAL))) +
  geom_density(fill = "purple", alpha = 1) +
  labs(
    x = "QUAL",
    y = "Density",
    title = "Distribution of QUAL"
  ) 
annotation_table[annotation_table$QUAL>50000000,]
ggplot(annotation_table, aes(x = QUAL)) + xlim(0,500) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 30, colour = "red") +
  gglayer_theme
  #geom_vline(xintercept = 2, colour = "blue")
cutoff <- quantile(annotation_table$QUAL, 0.05, na.rm = TRUE)


gglayer_theme <- list(
  theme_bw(),
  scale_color_brewer(palette="Dark2")
)

##QD 
#blue int is standard, red is actual five percent
QD_hist <- ggplot(annotation_table, aes(x = QD)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 9, colour = "red")+
  geom_vline(xintercept = 2, colour = "blue") +
  gglayer_theme
cutoff <- quantile(annotation_table$QD, 0.05, na.rm = TRUE)

##FS
FS_hist <- ggplot(annotation_table, aes(x = FS)) + xlim(0,100) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 18.5, colour = "red") +
  geom_vline(xintercept = 60, colour = "blue") +
  gglayer_theme
cutoff <- quantile(annotation_table$FS, 0.95, na.rm = TRUE)
#good at 5
##SOR
SOR_hist <- ggplot(annotation_table, aes(x = SOR)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 2.83, colour = "red") +
  geom_vline(xintercept = 3, colour = "blue") +
  gglayer_theme
cutoff <- quantile(annotation_table$SOR, 0.95, na.rm = TRUE)

##MQ
MQ_hist <- ggplot(annotation_table, aes(x = MQ)) + xlim(33,60) + ylim(0,10000) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 48.5, colour = "red") +
  geom_vline(xintercept = 40, colour = "blue") +
  gglayer_theme
cutoff <- quantile(annotation_table$MQ, 0.05, na.rm = TRUE)

##MQRankSum
MQRS_hist <- ggplot(annotation_table, aes(x = MQRankSum)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = -1.15, colour = "red") +
  geom_vline(xintercept = -12.5, colour = "blue") +
  gglayer_theme

cutoff <- quantile(annotation_table$MQRankSum, 0.05, na.rm = TRUE)

##ReadPosRankSum
RP <- ggplot(annotation_table, aes(x = ReadPosRankSum)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = -2.55, colour = "red") +
  geom_vline(xintercept = -8.0, colour = "blue") +
  gglayer_theme
cutoff <- quantile(annotation_table$ReadPosRankSum, 0.05, na.rm = TRUE)

hist(missing_indiv$F_MISS,
     breaks=50,
     col="blue",
     main="Sample Missingness Distribution",
     xlab="F_MISS (Proportion Missing)")

plot(density(missing_indiv$F_MISS),
     main="Density of Sample Missingness",
     xlab="F_MISS",
     lwd=2)

abline(v=0.1, col="red", lwd=2, lty=2)

sorted <- sort(missing_indiv$F_MISS)

plot(sorted,
     type="l",
     main="Sorted Sample Missingness",
     xlab="Sample rank",
     ylab="F_MISS")

abline(h=0.1, col="red", lty=2)

quantile(missing_indiv$F_MISS, c(0.9, 0.95, 0.99))
#       90%        95%        99% 
#0.0555556 0.0972222 0.4027780


hist(missing_site$F_MISS,
     breaks=50,
     col="blue",
     main="Sample Missingness Distribution",
     xlab="F_MISS (Proportion Missing)")

plot(density(missing_site$F_MISS),
     main="Density of Sample Missingness",
     xlab="F_MISS",
     lwd=2)

abline(v=0.1, col="red", lwd=2, lty=2)

sorted <- sort(missing_site$F_MISS)

plot(sorted,
     type="l",
     main="Sorted Sample Missingness",
     xlab="Sample rank",
     ylab="F_MISS")

abline(h=0.1, col="red", lty=2)

quantile(missing_site$F_MISS, c(0.9, 0.95, 0.99))
#90%       95%       99%  
#0.0555556 0.0972222 0.4027780

sum(missing_indiv$F_MISS > 0.05)
missing_indiv[missing_indiv$F_MISS > 0.05, ]
missing_site[missing_site$F_MISS > 0.05, ]
sum(missing_site$F_MISS > 0.05)

##calculate the depth per sample 

depth_table <- read.table("/Users/amritatrehan/documents/depth_table.tsv", header=TRUE)

# remove metadata columns
depth <- depth_table[, 4:ncol(depth_table)]

sample_means <- colMeans(depth, na.rm=TRUE)

sample_means

sample_means$Sample <- sub("\\.DP$", "", sample_means$Sample)

install.packages("devtools")
devtools::install_github("thomasp85/patchwork")

library(patchwork)
#Combination plot of the distribution of the quality score histogram plots, red is five percent tail that was cut off and blue is the standard
p1 <- FS_hist + MQ_hist + RP + MQRS_hist +SOR_hist + QD_hist + plot_annotation(title = "Distribution of Quality Scores", tag_levels = "A")
