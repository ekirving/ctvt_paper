#!/usr/bin/env Rscript
suppressWarnings(library("ape"))

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
data_file=args[1]
outgroup=args[2]
tree_file=args[3]
pdf_file=args[4]

# TODO remove when done testing
setwd("/Users/Evan/Dropbox/Code/ctvt")
data_file <- "njtree/test-pops.merged_map.geno.data"
outgroup <- "AndeanFox"
tree_file <- "njtree/test-pops.merged_map.geno.tree"
pdf_file <- "pdf/test-pops.merged_map.njtree.pdf"

m<-as.matrix(read.table(data_file, head=T, row.names=1))
tr=bionj(m)
tr<-root(tr, outgroup = outgroup, resolve.root = TRUE)

# sort the tree
tr <- ladderize(tr)

collist=c("#524FA1", "#FDB913", "lightgreen", "#00ADDC", "Darkgreen", "#ED1C24",
          "Black", "Pink", "Brown", "Cyan", "midnightblue", "palevioletred3",
          "lightcoral", "yellow4", "wheat4")

# TODO colour the tips of the tree
# col <- gsub("[^A-Z]", "", substring(tr$tip, 0, 3))
# col[col=="DOM"] = collist[1]
# col[col=="FRE"] = collist[2]
# col[col=="IB1"] = collist[3]
# col[col=="IB2"] = collist[4]
# col[col=="OUT"] = collist[5]

# save the tree data
write.tree(tr, file=tree_file)

# calcualte the plot size
numnodes <- length(m[,1])
plotsize = max(c(numnodes/20, 7))

# plot the tree
# pdf(file=pdf_file, width = plotsize, height = plotsize)
plot(tr, type='fan', cex=0.8) #, tip.color=col)
# dev.off()
