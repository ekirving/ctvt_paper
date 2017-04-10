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
# pca1_file = "smartpca/all-pops.ascertain1.calc.pca"
# pca2_file = "smartpca/all-pops.ascertain1.proj.pca"
# pve_file = "smartpca/all-pops.ascertain1.pve"
# pdf_file = "pdf/all-pops.ascertain1.PCA.1.2.pdf"
# comp1=1
# comp2=2
# labeled=strtoi("0")

# get the percentage of variance each component explains
pve <- round(read.table(pve_file)[,1]*100, 1)

# read in the PCA data
t1 <- read.table(pca1_file)
t2 <- read.table(pca2_file)

# join the two data frames
dat <- rbind(t1, t2)
names(dat)[1] <- "Code"

# get the metadata matrix for all the samples
info <- read.table('pop_names.csv', sep = ",", quote = '"', header=TRUE)

# join the population types
total <- merge(dat, info[c(1,3)], by = 'Code')

# setup the different symbol types
# see http://www.cookbook-r.com/Graphs/Shapes_and_line_types/ for shape codes
shapes <- c(15, 17, 19, 18, 13, 9, 0)

# make sure there are enough of each shape
shapes <- rep_len(shapes, length.out=length(unique(total$Type.Name)))

# change the order of the factors
# total[,'Type.Name'] <- factor(total[,'Type.Name'], levels = c("Taurine, Eurasian", "Taurine, African", "Zebu", "Hybrid", "Outgroup", "Auroch", "Socotra"))

# setup the colours
colours <- c('#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02')
colours <- rep_len(colours, length.out=length(unique(total$Type.Name)))

alpha=c(1, 1, 1, 0.1, 1, 1, 1, 1, 1, 1, 1, 1)
pdf(file=pdf_file, width = 10, height = 7)

gg <- ggplot(total, aes(total[[comp1+2]], total[[comp2+2]])) +
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
  gg <- gg + geom_text(aes(label=sub("^[^:]*:", "", dat$V2)), hjust=-.3, vjust=0)
}


# display the plot
gg

dev.off()
