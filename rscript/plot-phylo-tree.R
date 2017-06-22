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
# data_file <- "njtree/all-pops.merged_v1.geno.data"
# treetype <- "phylogram" # "phylogram" / "fan"
# outgroup <- "AndeanFox"
# tree_file <- "njtree/all-pops.merged_v1.geno.tree"
# pdf_file <- "pdf/all-pops.merged_v1.geno.njtree.pdf"

# load the distance matrix
m<-as.matrix(read.table(data_file, head=T, check.names=FALSE))

# build the NJ tree
tr <- bionj(m[,-c(1:2)])

# root the tree
tr <- root(tr, outgroup = outgroup, resolve.root = TRUE)

# resolve.root adds a zero-length branch below the MRCA of the ingroup, so lets
# find the outgroup branch and place the root in the middle of the branch
last <- length(tr$edge.length)
len <- tr$edge.length[last]

tr$edge.length[1] <- len/2     # the edge to the MRCA is always first
tr$edge.length[last] <- len/2  # the edge to the outgroup is always last

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
meta <- merge(m, info[c(1,3,4)], by = 'Code')
rownames(meta) <- meta[,'Sample']

# get the population names and colours for the legend
key<-unique(meta[c('Type.Name','Colour')])
# key<-key[with(key, order(Type.Name)), ]
key[] <- lapply(key, as.character)

# fix the type issue with the colour codes
cols <- sapply(meta[tr$tip.label, 'Colour'], as.character)

pdf(file=pdf_file, width = plotsize, height = plotsize)

# plot the tree
plot(tr, type=treetype, cex=0.6, tip.color=cols, label.offset=0.001)

# add the legend
legend("bottomleft", legend = key[[1]], fill = key[[2]], cex=0.6)

dev.off()
