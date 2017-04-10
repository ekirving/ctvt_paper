#!/usr/bin/env Rscript
suppressWarnings(library(ggplot2))
suppressWarnings(library(scales))
library(stringr)

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
pca1_file=args[1]
pca2_file=args[2]
pve_file=args[3]
pdf_file=args[4]
comp1=strtoi(args[5])
comp2=strtoi(args[6])
labeled=strtoi(args[7])

# TODO remove when done testing
# setwd("/Users/Evan/Dropbox/Code/socotra")
# pca1_file = "smartpca/all-pops.calc.pca"
# pca2_file = "smartpca/all-pops.proj.pca"
# pve_file = "smartpca/all-pops.pve"
# pdf_file = "pdf/all-pops.PCA.1.2.pdf"
# comp1=1
# comp2=2
# labeled=strtoi("0")

# get the percentage of variance each component explains
pve <- round(read.table(pve_file)[,1]*100, 1)

# read in the PCA data
t1 = read.table(pca1_file)
names(t1)[1]<-paste("Population")

t2 = read.table(pca2_file)
names(t2)[1]<-paste("Population")

# setup the shapes, see http://www.cookbook-r.com/Graphs/Shapes_and_line_types/ for shape codes
solid=c(15, 17, 15, 20, 18, 17, 15, 16, 20, 18, 15) # solid shapes
hollow=c(21, 22, 23, 24, 25) # hollow shapes (projected)

# make sure there are enough of each shape
v=c(rep_len(solid, length.out=length(unique(t1$Population))), rep_len(hollow, length.out=length(unique(t2$Population))))

# join the two data frames
dat=rbind(t1, t2)

# setup the colours
col=c("#524FA1", "#FDB913", "lightgreen", "#00ADDC", "Darkgreen", "#ED1C24", "Black", "Brown", "lightcoral", "midnightblue", "palevioletred3", "Cyan", "yellow4", "wheat4")
col=rep_len(col, length.out=length(unique(dat$Population)))

alpha=c(1, 1, 1, 0.1, 1, 1, 1, 1, 1, 1, 1, 1)
pdf(file=pdf_file, width = 10, height = 7)
gg <- ggplot(dat, aes(dat[[paste('V', comp1+2, sep='')]], dat[[paste('V', comp2+2, sep='')]])) + aes(shape=factor(Population)) + scale_shape_manual(values=v) + geom_point(aes(colour = factor(Population)), size=3, alpha=1) +  xlab(paste("PC", comp1, " (", pve[comp1], "%)", sep='')) + ylab(paste("PC", comp2, " (", pve[comp2], "%)", sep='')) + scale_colour_manual(values=col) + theme_bw() + theme(legend.title=element_blank(), legend.key = element_blank()) + guides(colour = guide_legend(override.aes = list(size=4)))
if (labeled) {
  # label all the points
  gg <- gg + geom_text(aes(label=str_sub(dat$V2, -3)),hjust=-.3, vjust=0)
}
# display the plot
gg
dev.off()