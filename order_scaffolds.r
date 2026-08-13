####FEMALE MAPS #### 

##CHROM 1##

#Get ouput from ordermarkers2 after snp.txt

amrita_folder <- "/Users/amritatrehan/documents/Thesis/Marey_maps/no_small_5/"
james_folder <- "/Users/jamescotton/Work_newlaptop/Glasgow_Teaching/BIOINF_MSc/Amrita_2026/Map_ordering"
#you'll need to change this
setwd(amrita_folder)
#and swap this back
#fai_name <- "/Users/jamescotton/Work_newlaptop/Glasgow_Teaching/BIOINF_MSc/Amrita_2026/Map_ordering/dracunculus_medinensis_v3.fa.fai"
fai_name <- "/Users/amritatrehan/documents/Thesis/dracunculus_medinensis_v3.fa.fai"

#can you write a loop through chromosomes? 
#
#if the filenames are "simple you can just do:
#for (chrom_num in 1:3) { 
# then use paste0 to generate the map_file_name and the output file names
#}
##  
  #map_file_name <- paste0(
  #  "order", chrom_num, "_eval8_order.mapped"
  #)
map_file_name = "order1_eval8_order.mapped"
#map_file_name = "order1_eval_male_11_order.mapped.fixed"

#the stuff below has been rewritten to use a column specified by this string
map_column = "female_cM"

#--- try to write some re-usable code --

for (chrom_num in 1:3) { 
  map_file_name <- paste0(
    "order", chrom_num, "_eval8_order.mapped"
  )
  map1 <- read.table(
    map_file_name,
    header = FALSE,
    comment.char = "#",
    stringsAsFactors = FALSE
  )
  
  names(map1)[1:4] <- c("scaffold", "bp", "male_cM", "female_cM")
  #needed for some reason..
  map1$bp <- as.numeric(map1$bp)
  
  #male or female map depending on value of map_column variable
  
  library(ggplot2)
  library(dplyr)
  library(stringr)
  #use a string variable to refer to column name..
  ggplot(map1, aes(x = bp, y = .data[[map_column]])) +
    geom_point(size = 1) +
    facet_wrap(~ scaffold, scales = "free_x") +
    labs(
      x = "Physical position (bp)",
      y = paste0(str_replace(map_column,"_cM","")," genetic position (cM)")
    ) +
    theme_bw()
  #make filename from variable 
  ggsave(paste0(str_replace(map_column,"_cM",""),"_facet_map_chrom1.pdf"), width = 18, height = 14)
  
  #get fai file with contig lenghts (scaf_lnegths)
  fai <- read.table(fai_name, sep = "\t")
  names(fai)[1:3] <- c("scaffold", "length", "size")
  
  #create a table with scaf_lengths, their names, and the number of markers (SNPS)
  #re-written to pass column as string
  scaf_summary <- map1 %>%
    left_join(fai[, c("scaffold", "length")], by = "scaffold") %>%
    group_by(scaffold) %>%
    summarise(
      length = first(length),
      mean_map = mean(!!as.symbol(map_column), na.rm = TRUE),
      n_markers = n()
    ) %>%
    arrange(mean_map)
  
  # rearrange scaf_summary so it is in order of mean_map 
  #alreay in that order..
  #scaf_summary
  #scaf_summary <- scaf_summary %>%
  #  arrange(mean_map)
  
  # calculate start using cum length minus this contigs length to get to the start, amd 0 to begin 
  
  scaf_summary <- scaf_summary %>%
    mutate(
      start = c(0, cumsum(length)[-length(length)])
    )
  
  #make ordered graph (no flips yet) and add scaffold starts 
  
  ordered_map <- map1 %>%
    left_join(
      scaf_summary %>% select(scaffold, start),
      by = "scaffold"
    )
  
  ordered_map <- ordered_map %>%
    mutate(fake_bp = start + bp)
  
  #find which sections need to be flipped 
  
  orientation <- map1 %>%
    group_by(scaffold) %>%
    summarise(
      n = n(),
      cor = if (n() < 2 ||
                sd(bp) == 0 ||
                sd(!!as.symbol(map_column)) == 0) {
        NA_real_
      } else {
        cor(bp, !!as.symbol(map_column))
      }
    ) %>%
    mutate(
      flip = cor < 0 & n >= 5
    )
  
  #which should be flipped 
  
  orientation %>%
    filter(flip)
  
  #which should be looked at 
  
  orientation %>%
    filter(is.na(cor))
  
  #calc new starting
  #flip or not and add to map
  
  map_flip <- map1 %>%
    left_join(orientation, by = "scaffold") %>%
    left_join(fai[, c("scaffold", "length")], by = "scaffold") %>%
    mutate(
      bp_flip = if_else(flip, length - bp, bp)
    )
  
  fai2 <- fai %>%
    rename(scaf_length = length)
  
  #re-do scaf_summary with flipped contigs
  
  scaf_summary <- map_flip %>%
    group_by(scaffold) %>%
    summarise(
      scaf_length = first(length),
      mean_map = mean(!!as.symbol(map_column), na.rm = TRUE),
      n_markers = n()
    ) %>%
    arrange(mean_map)
  
  scaf_summary <- scaf_summary %>%
    mutate(offset = c(0, cumsum(head(scaf_length, -1))))
  
  #order them and create fake genomic coordinates for flipped
  
  map_ordered <- map_flip %>%
    left_join(scaf_summary[, c("scaffold", "offset")], by = "scaffold") %>%
    mutate(fake_bp = bp_flip + offset)
  
  flip_map <- ggplot(map_ordered, aes(fake_bp, !!as.symbol(map_column), colour = scaffold)) +
    geom_point(size = 1) +
    theme_bw()
  
  assign(
    paste0(str_replace(map_column, "_cM", ""), "map", chrom_num),
    flip_map
  )
  
  ggsave(paste0(str_replace(map_column,"_cM",""),"_facet_map_chrom", chrom_num, ".flipped.pdf"), plot=flip_map, width = 18, height = 14)
  
  ####MAREY MAP ####
  library(tidyr)

  #getting confidence intervals for chromosome-wide recombination rate estimtes
  model_formula1 = paste0("fake_bp ~ ",map_column)
  print(model_formula1)
  fit <- lm(data=map_ordered,as.formula(model_formula1))
  #The 'Estimate' falue here for male_cM /female_cM is bp per cM
  summary(fit)
  confint(fit)
  model_formula2 = paste0(map_column," ~ I(fake_bp/1000000)")
  print(model_formula2)
  fit <- lm(data=map_ordered,as.formula(model_formula2))
  #The 'Estimate' falue here for fake_bp is not in cM per Mb
  summary(fit)
  confint(fit)

  #get chrom length, # of markers, and genetic map length in cM, Recombination rate in males and females   
  chrom_length <- max(map_ordered$fake_bp, na.rm = TRUE)
  chrom_length
  n_markers <- sum(!is.na(map_ordered$female_cM))
  n_markers
  map_length_cM <- max(map_ordered$female_cM, na.rm = TRUE) -
    min(map_ordered$female_cM, na.rm = TRUE)
  map_length_cM
  #doing LOESS for Marey Map
  
  #span here will set the smoothness of the loess
  #you might want to use a lower number than this, but then its more obvious that recombination rate
  #goes negative at the strange bump in your map. A bigger number makes the loess fit worse but disguises the -ve a bit.
  wiggle <- loess(data=map_ordered,model_formula2,span=0.5,degree=2)
  #the 1000,000 here is the number of points used in ploting the smoother, 
  #but also the size of windows used in recombinatuion rate regression in the next bit
  fake_data <- data.frame(fake_bp=seq(100000,max(na.omit(map_ordered$fake_bp)),by=100000))
  fake_data$cM <- predict(wiggle,newdata=fake_data)
  
  #plot smoother
  ggplot(data=map_ordered,aes(x=fake_bp,y=!!as.symbol(map_column), colour = scaffold)) +
    geom_point(size = 1) + geom_line(data=fake_data,mapping=aes(x=fake_bp,y=cM),inherit.aes = FALSE) + 
    theme_bw()
  
  #make recombination rate
  #use slider package, and fit to window of 10 values in loess
  #could also do this direclty on the data (not the loess) but would be a bit more fiddly to code!
  
  #for slide function
  library(slider)
  #for tidy
  library(broom)
  #for map
  library(purrr)
  
  #before and after here set the size of the window over which
  #recombination rate is estimated - in terms of the number of windows of size set in fake_bp
  
  fake_data_plus_recomb <- fake_data %>% mutate(Mb = fake_bp/1E6 )
  fake_data_plus_recomb$models <- slide(
    fake_data_plus_recomb, 
    ~lm(cM ~ Mb, data = .x), 
    .before = 10,
    .after = 10,
    .complete = TRUE
  )
  fake_data_plus_recomb <- fake_data_plus_recomb %>% mutate(tidied_models = purrr::map(models,tidy)) %>% unnest(tidied_models)
  fake_data_plus_recomb <- fake_data_plus_recomb %>% filter(term == "Mb")
  
  
  
  #plot smoother AND recombination rate 
  #calculations to fool ggplot2 into plotting two different y axes
  ylim.prim <- range(map_ordered[[map_column]])
  ylim.sec <- range(fake_data_plus_recomb$estimate)
  b <- diff(ylim.prim)/diff(ylim.sec)
  a <- ylim.prim[1] - b*ylim.sec[1]
  graph <- ggplot(data=map_ordered,aes(x=fake_bp,y=!!as.symbol(map_column), colour = scaffold)) +
    geom_point(size = 1) + geom_line(data=fake_data_plus_recomb,mapping=aes(x=fake_bp,y=cM),inherit.aes = FALSE) + 
    geom_line(data=fake_data_plus_recomb,mapping=aes(x=fake_bp,y=estimate*10), inherit.aes=FALSE, color="red") + 
    scale_y_continuous(sec.axis = sec_axis(~ (. - a)/b , name = "recombination rate cM/Mb")) + 
    theme_bw() + theme( axis.title.y.right = element_text(color="red") ,axis.text.y.right = element_text(color = "red")  ) 
  
  assign(
    paste0(str_replace(map_column, "_cM", ""), "graph", chrom_num),
    graph
  )
  
  #loess.ci
  library(spatialEco)
  map_ordered.complete <- na.omit(map_ordered)
  wiggle.ci <- loess.ci(y=map_ordered.complete[[map_column]],x=map_ordered.complete$fake_bp,span=0.5)
  #the 1000,000 here is the number of points used in ploting the smoother, 
  #but also the size of windows used in recombinatuion rate regression in the next bit
  plot_frame <- data.frame(bp=map_ordered.complete$fake_bp,scaffold=map_ordered.complete$scaffold,cM=map_ordered.complete[[map_column]],smooth_cM=wiggle.ci$loess,lci=wiggle.ci$lci,uci=wiggle.ci$uci)
  loess_plot <- ggplot(plot_frame,aes(x=bp,y=cM,color=scaffold)) + geom_point() + geom_line(aes(y=smooth_cM),color="black") + geom_ribbon(aes(ymin=lci,ymax=uci),color="grey",alpha=0.5) + theme_bw()
  
  assign(
    paste0(str_replace(map_column, "_cM", ""), "loess_plot", chrom_num),
    loess_plot
  )
}

library(patchwork)
p1 <- femalemap1 + ggtitle ("female chromosome 1")
p2 <- femalemap2 + ggtitle ("female chromosome 2")
p3 <- femalemap3 + ggtitle ("female chromosome 3")
p4 <- malemap1 + ggtitle ("male chromosome 1")
p5 <- malemap2 + ggtitle ("male chromosome 2")
p6 <- malemap3 + ggtitle ("male chromosome 3")

ggsave(("female_vs_male_maps.pdf"), plot=(p1 | p4 ) /
         (p2 | p5) /
         (p3 | p6), width = 18, height = 14)

ggsave(("female_maps.pdf"), plot=(p1 | p2 ) /
         (p3 | plot_spacer()), width = 18, height = 14)

ggsave(("male_maps.pdf"), plot=(p4 | p5 ) /
         (p6 | plot_spacer()), width = 18, height = 14)

##maray map plots 
p7 <- femalegraph1 + ggtitle ("female chromosome 1") + theme(legend.position = "none")
p8 <- femalegraph2 + ggtitle ("female chromosome 2") + theme(legend.position = "none")
p9 <- femalegraph3 + ggtitle ("female chromosome 3") + theme(legend.position = "none")
p10 <- femaleloess_plot1 + ggtitle ("female chromosome 1") + theme(legend.position = "none")
p11 <- femaleloess_plot2 + ggtitle ("female chromosome 2") + theme(legend.position = "none")
p12 <- femaleloess_plot3 + ggtitle ("female chromosome 3") + theme(legend.position = "none")

#female maps

ggsave(("female_maray.pdf"), plot=(p7 | p8 ) /
         (p9 | plot_spacer()), width = 18, height = 14)

ggsave(("female_loess.pdf"), plot=(p10 | p11 ) /
         (p12 | plot_spacer()), width = 18, height = 14)

#male maps

p13 <- malegraph1 + ggtitle ("male chromosome 1") + theme(legend.position = "none")
p14 <- malegraph2 + ggtitle ("male chromosome 2") + theme(legend.position = "none")
p15 <- malegraph3 + ggtitle ("male chromosome 3") + theme(legend.position = "none")
p16 <- maleloess_plot1 + ggtitle ("male chromosome 1") + theme(legend.position = "none")
p17 <- maleloess_plot2 + ggtitle ("male chromosome 2") + theme(legend.position = "none")
p18 <- maleloess_plot3 + ggtitle ("male chromosome 3") + theme(legend.position = "none")

#female maps

ggsave(("male_maray.pdf"), plot=(p13 | p14 ) /
         (p15 | plot_spacer()), width = 18, height = 14)

ggsave(("male_loess.pdf"), plot=(p16 | p17 ) /
         (p18 | plot_spacer()), width = 18, height = 14)

####average rate per chrom + confidence ####

male_fit1 <- lm(malemap1$data$male_cM ~ malemap1$data$fake_bp)
#get estimate
summary(male_fit1)
#get confidence interval
confint(male_fit1)
