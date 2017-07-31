#!/usr/bin/env Rscript
library("ape")

# get the command line arguments
# args <- commandArgs(trailingOnly = TRUE)
# data_file=args[1]
# treetype=args[2]
# outgroup=args[3]
# tree_file=args[4]
# pdf_file=args[5]

# TODO remove when done testing
setwd("/Users/Evan/Dropbox/Code/ctvt")
data_file <- "njtree/all-pops.merged_v2.geno.data"
treetype <- "phylogram" # "phylogram" / "fan"
outgroup <- "AndeanFox"

tree_file <- '/Users/Evan/Downloads/trees/merged_map_tv_standard.nex.con.tre'
pdf_file  <- '/Users/Evan/Downloads/trees/pdf/merged_map_tv_standard.nex.con.pdf'

# tree_file <- '/Users/Evan/Downloads/trees/merged_map_tv_standard1.nex.con.tre'
# pdf_file  <- '/Users/Evan/Downloads/trees/pdf/merged_map_tv_standard1.nex.con.pdf'

# tree_file <- '/Users/Evan/Downloads/trees/merged_pruned_bootstrap_VIET.tree'
# pdf_file  <- '/Users/Evan/Downloads/trees/pdf/merged_pruned_bootstrap_VIET.pdf'

# tree_file <- '/Users/Evan/Downloads/trees/merged_pruned_bootstrap.tree'
# pdf_file  <- '/Users/Evan/Downloads/trees/pdf/merged_pruned_bootstrap.pdf'

# load the distance matrix
m<-as.matrix(read.table(data_file, head=T, check.names=FALSE))

# load the NJ tree from file
# tr <- read.tree(tree_file)
tr <- read.nexus(tree_file)

# root the tree
tr <- root(tr, outgroup = outgroup, resolve.root = TRUE, edgelabel = TRUE)

# fix the node labels
if (!is.null(tr$node.label)) {
    tr$node.label <- strtrim(tr$node.label, 4)
    tr$node.label[1] <- ''
}

# resolve.root adds a zero-length branch below the MRCA of the ingroup, so lets
# find the outgroup branch and place the root in the middle of the branch
last <- length(tr$edge.length)
len <- tr$edge.length[last]

tr$edge.length[1] <- len/2     # the edge to the MRCA is always first
tr$edge.length[last] <- len/2  # the edge to the outgroup is always last

# sort the tree
tr <- ladderize(tr)

# calcualte the plot size
numnodes <- length(m[,1])
plotsize = max(c(numnodes/20, 7))

# get the metadata matrix for all the samples
info <- read.table('pop_names.csv', sep = ",", quote = '"', header=TRUE, comment.char="")

# join the metadata, so we can get the node colours
meta <- merge(m, info[c('Code','Type.Name','Colour','Order')], by = 'Code')
rownames(meta) <- meta[,'Sample']

# get the population names and colours for the legend
key<-unique(meta[c('Type.Name','Colour','Order')])
key<-key[with(key, order(Order)), ]
key[] <- lapply(key, as.character)

# fix the type issue with the colour codes
cols <- sapply(meta[tr$tip.label, 'Colour'], as.character)

pdf(file=pdf_file, width = plotsize, height = plotsize)

# plot the tree
if (is.null(tr$node.label)) {
    plot(tr, type=treetype, cex=0.6, tip.color=cols, label.offset=0.001, show.node.label=F)
} else {
    plot(tr, type=treetype, cex=0.6, tip.color=cols, label.offset=0.001, show.node.label=T)
}
legend(x=0,y=25,legend = key[[1]], fill = key[[2]], cex=0.6)
add.scale.bar(cex = 0.6)

dev.off()
