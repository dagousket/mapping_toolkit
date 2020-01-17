#!/usr/bin/env Rscript
# This script makes an overlap summary file for three coordinate bed files.
# Help with Rscript --vanilla summary_region_overlaps.R --help

library("optparse")
library("tidyverse")
library("GenomicRanges")
library("UpSetR")
library("eulerr")
library("RColorBrewer")
library("gridExtra")
library("grid")

option_list = list(
  make_option(c("--bed1"), type="character",  
              help="First BED file", metavar="character"),
  make_option(c("--bed2"), type="character",  
              help="Second BED file", metavar="character"),
  make_option(c("--bed3"), type="character",  
              help="Third BED file", metavar="character"),
  make_option(c("--out"), type="character", 
              help="path and name of the output PDF", default = './results.pdf')
); 

opt_parser = OptionParser(option_list=option_list, description = "An R script to perfrom overlaps between 3 BED files and extract main features.");
opt = parse_args(opt_parser);


if (is.null(opt$bed1)){
  print_help(opt_parser)
  stop("Some arguments are missing", call.=FALSE)
}


## Make overlap stats based on 3 sets of genomic regions and eventually an extra gene one.

# Example of how intersection are build :
#
# BED1    : - - ======= - - - - ==============- - - - - - - - - - ==- -
# BED2    : - ======- - - - - - -=====- - - === - - - - - - - - - - - -
# BED3    : - - - ==================- - - - - - - - - ====- - - - - - -
#
# Overlap : - - - ==- - - - - - -===- - - - ==- - - - ====- - - - ==- -
# Depth   :       3               3          2          1          1

bed1 <- read.table(file = opt$bed1)
bed2 <- read.table(file = opt$bed2)
bed3 <- read.table(file = opt$bed3)

gr1 <- GRanges(seqnames = Rle(bed1$V1), ranges = IRanges(start = bed1$V2, end = bed1$V3, names = bed1$V4))
gr2 <- GRanges(seqnames = Rle(bed2$V1), ranges = IRanges(start = bed2$V2, end = bed2$V3, names = bed2$V4))
gr3 <- GRanges(seqnames = Rle(bed3$V1), ranges = IRanges(start = bed3$V2, end = bed3$V3, names = bed3$V4))


# Build set of 2 layers and 3 layers regions
gr_1.2 <- setdiff(intersect(gr1, gr2), gr3)
mcols(gr_1.2)$score <- 2
mcols(gr_1.2)$overlap <- "1.2"
gr_2.3 <- setdiff(intersect(gr2, gr3), gr1)
mcols(gr_2.3)$score <- 2
mcols(gr_2.3)$overlap <- "2.3"
gr_1.3 <- setdiff(intersect(gr1, gr3), gr2)
mcols(gr_1.3)$score <- 2
mcols(gr_1.3)$overlap <- "1.3"
gr_1.2.3 <- intersect(intersect(gr1, gr2),gr3)
mcols(gr_1.2.3)$score <- 3
mcols(gr_1.2.3)$overlap <- "1.2.3"

# Build set of non-overlapping regions
gr_1alone <- setdiff(setdiff(gr1, gr2), gr3)
mcols(gr_1alone)$score <- 1
mcols(gr_1alone)$overlap <- "1"
gr_2alone <- setdiff(setdiff(gr2, gr3), gr1)
mcols(gr_2alone)$score <- 1
mcols(gr_2alone)$overlap <- "2"
gr_3alone <- setdiff(setdiff(gr3, gr1), gr2)
mcols(gr_3alone)$score <- 1
mcols(gr_3alone)$overlap <- "3"

# Make a function to expand compress bedgraph format to single base pair format
expand_bed <- function(df){
  df_out <- c()
  for (i in 1:nrow(df)){
    df_out <- rbind(df_out, data.frame('seqname' = df[i,'seqnames'], 'bp' = seq(df[i,'start'],df[i,'end']),'score' = df[i,'score'],'overlap'=df[i, 'overlap']))
  }
  return(df_out)
}

# Combine all sets
gr <- as.data.frame(c(gr_1alone, gr_2alone, gr_3alone, gr_1.2, gr_1.3, gr_2.3, gr_1.2.3))
gr <- gr[order(gr$start),]
gr_expanded <- expand_bed(gr)

# Plot frequency of overlap type and proportion of region overlap
gr_expanded <- gr_expanded %>% mutate('layer1' = ifelse(grepl('1',overlap), 1, 0), 'layer2' = ifelse(grepl('2',overlap), 1, 0), 'layer3' = ifelse(grepl('3',overlap), 1, 0))
upset(gr_expanded, sets = c("layer1","layer2","layer3"), sets.bar.color = "#56B4E9", order.by = "freq")
grid.edit('arrange',name='arrange3')
upset_plot = grid.grab()

layer_color = colorRampPalette(brewer.pal(4,"Dark2"))(4)[c(3,2,1,4)]
fit1 <- euler(gr_expanded[,c("layer1","layer2","layer3")])
euler_plot <- plot(fit1,
   fill = layer_color[1:3],
   border = "transparent",
   fill_opacity = 0.6,
   cex = 1,
   labels = c('bed1','bed2','bed3'),
   counts = FALSE)
  
pdf(file = opt$out, onefile = TRUE, width = 15 , height = 7)
lay <- matrix(c(1,2), nrow = 1)
grid.arrange(arrangeGrob(upset_plot, top = 'Upset graph'),arrangeGrob(euler_plot, top = 'Venn diagram'), layout_matrix = lay)
dev.off()
