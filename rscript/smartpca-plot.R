#!/usr/bin/env Rscript
suppressWarnings(library(ggplot2))
suppressWarnings(library(scales))
suppressWarnings(library(stringr))

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
# setwd("/Users/Evan/Dropbox/Code/ctvt")
# pca1_file = "smartpca/all-pops.merged_map.prj-DPC.calc.pca"
# pca2_file = "smartpca/all-pops.merged_map.prj-DPC.proj.pca"
# pve_file = "smartpca/all-pops.merged_map.prj-DPC.pve"
# pdf_file = "pdf/all-pops.merged_map.PCA.1.2.pdf"
# comp1=1
# comp2=2
# labeled=strtoi("0")

# get the percentage of variance each component explains
pve <- round(read.table(pve_file)[,1]*100, 1)

# get the metadata matrix for all the samples
info <- as.data.frame(read.table('pop_names.csv', sep = ",", quote = '"', header=TRUE, comment.char=""))

# read in the PCA data
t1 <- read.table(pca1_file)
t2 <- read.table(pca2_file)

names(t1)[1] <- "Code"
names(t2)[1] <- "Code"

# join the population types
t1meta <- merge(t1, info[c(1,3,4)], by = 'Code')
t2meta <- merge(t2, info[c(1,3,4)], by = 'Code')

# join the two data frames
meta <- rbind(t1meta, t2meta)

# setup the different symbol types
# see http://www.cookbook-r.com/Graphs/Shapes_and_line_types/ for shape codes
solid=c(15, 17, 15, 20, 18, 17, 15, 16, 20, 18, 15) # solid shapes
hollow=c(21, 22, 23, 24, 25) # hollow shapes (projected)

# make sure there are enough of each shape
shapes<-c(rep_len(solid, length.out=length(unique(t1meta$Type.Name))),
          rep_len(hollow, length.out=length(unique(t1meta$Type.Name))))

# setup the colours
colours <- unique(as.character(meta$Colour))

alpha=c(1, 1, 1, 0.1, 1, 1, 1, 1, 1, 1, 1, 1)
pdf(file=pdf_file, width = 10, height = 7)

gg <- ggplot(meta, aes(meta[[comp1+2]], meta[[comp2+2]])) +
    aes(shape=factor(Type.Name)) +
    scale_shape_manual(values=shapes) +
    geom_point(aes(colour = factor(Type.Name)), size=4, alpha=1) +
    xlab(paste("PC", comp1, " (", pve[comp1], "%)", sep='')) +
    ylab(paste("PC", comp2, " (", pve[comp2], "%)", sep='')) +
    scale_colour_manual(values=colours) +
    theme_bw() +
    # coord_fixed() +
    theme(legend.title=element_blank(), legend.key = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=4)))

if (labeled) {
  # label all the points
  gg <- gg + geom_text(aes(label=sub("^[^:]*:", "", meta$V2)), hjust=-.3, vjust=0, size=3)
}


# display the plot
gg

dev.off()
