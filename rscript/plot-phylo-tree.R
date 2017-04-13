#!/usr/bin/env Rscript
suppressWarnings(library("ape"))

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
data_file=args[1]
treetype=args[2]
outgroup=args[3]
tree_file=args[4]
pdf_file=args[5]

# TODO remove when done testing
# setwd("/Users/Evan/Dropbox/Code/ctvt")
# data_file <- "njtree/test-pops.merged_map.geno.data"
# treetype <- "phylogram" #"fan"
# outgroup <- "AndeanFox"
# tree_file <- "njtree/test-pops.merged_map.geno.tree"
# pdf_file <- "pdf/test-pops.merged_map.njtree.pdf"

m<-as.matrix(read.table(data_file, head=T, check.names=FALSE))
tr=bionj(m[,-c(1:2)])
tr<-root(tr, outgroup = outgroup, resolve.root = TRUE)

# sort the tree
tr <- ladderize(tr)

# save the tree data
write.tree(tr, file=tree_file)

# calcualte the plot size
numnodes <- length(m[,1])
plotsize = max(c(numnodes/20, 7))

# get the metadata matrix for all the samples
info <- read.table('pop_names.csv', sep = ",", quote = '"', header=TRUE, comment.char="")

# join the metadata, so we can get the node colours
meta <- merge(m, info[c(1,4)], by = 'Code')
rownames(meta) <- meta[,'Sample']

# fix the type issue with the colour codes
cols <- sapply(meta[tr$tip.label, 'Colour'], as.character)

# plot the tree
pdf(file=pdf_file, width = plotsize, height = plotsize)
plot(tr, type=treetype, cex=0.6, tip.color=cols)
# plot(tr, type=treetype, cex=0.8, tip.color=collist)
dev.off()
