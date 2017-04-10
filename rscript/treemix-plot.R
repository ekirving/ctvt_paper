#!/usr/bin/env Rscript
source("rscript/treemix-plotting_funcs.R")

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
data_file=args[1]
poplist_file=args[2]
pdf_file=args[3]

# TODO remove when done testing
# setwd("/Users/Evan/Dropbox/Code/socotra")
# data_file <- "treemix/select-pops.BovineHD.geno.k1000.m3"
# poplist_file <- "treemix/select-pops.BovineHD.poplist"
# pdf_file <- "pdf/select-pops.BovineHD.treemix.k1000.m3.pdf"

pdf(file=pdf_file, width = 16, height = 8)
par(mfrow=c(1,2))
plot_tree(data_file, poplist_file)
plot_resid(data_file, poplist_file)
dev.off()
