#!/usr/bin/env Rscript
suppressWarnings(library(ggplot2))

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
data_file <- args[1]
pdf_file <- args[2]
xlab <- args[3]
ylab <- args[4]

# TODO remove when done testing
# data_file  <-  "admix/select-pops.BovineSNP50v1.pruned.CV.data"
# pdf_file  <-  "pdf/select-pops.BovineSNP50v1.admix.CV.pdf"
# xlab <- "Ancestral populations (K)"
# ylab <- "Cross-validation Error"

# read the data file
dat <- read.table(data_file, header=TRUE)

# TODO highlight the 3 lowest points

pdf(file=pdf_file, width=10, height=7)

ggplot(dat, aes(x=dat[[1]], y=dat[[2]], group=1)) +
    geom_point() +
    geom_line(stat='identity') +
    theme(legend.title=element_blank(),
          legend.key = element_blank(),
          panel.background = element_blank()) +
    xlab(xlab) +
    ylab(ylab) +
    scale_x_continuous(breaks=c(0:nrow(dat))) +
    scale_y_continuous(breaks=seq(0, 3, 0.1))

dev.off()
