#!/usr/bin/env Rscript
suppressWarnings(library("ape"))

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
data_file=args[1]
outgroup=args[2]
tree_file=args[3]
pdf_file=args[4]

# TODO remove when done testing
# setwd("/Users/Evan/Dropbox/Code/ctvt")
# data_file <- "njtree/test-pops.merged_map.geno.data"
# outgroup <- "AndeanFox"
# tree_file <- "njtree/test-pops.merged_map.geno.tree"
# pdf_file <- "pdf/test-pops.merged_map.njtree.pdf"

m<-as.matrix(read.table(data_file, head=T, row.names=1))
tr=bionj(m)
tr<-root(tr, outgroup = outgroup, resolve.root = TRUE)

# sort the tree
tr <- ladderize(tr)

# # grab the tip names
# tips <- as.data.frame(tr$tip)
# names(tips)[1] <- 'Samples'
#
# # get the metadata matrix for all the samples
# info <- read.table('pop_names.csv', sep = ",", quote = '"', header=TRUE, comment.char="")
#
# # join the metadata
# total <- merge(dat, info[c(1,3)], by = 'Code')
#
# # setup the different symbol types
# # see http://www.cookbook-r.com/Graphs/Shapes_and_line_types/ for shape codes
# shapes <- c(15, 17, 19, 18, 13, 9, 0)
#
# # make sure there are enough of each shape
# shapes <- rep_len(shapes, length.out=length(unique(total$Type.Name)))
#
# # change the order of the factors
# # total[,'Type.Name'] <- factor(total[,'Type.Name'], levels = unique(total[,'Type.Name']))
#
# # setup the colours
# colours <- c('#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02')
# colours <- rep_len(colours, length.out=length(unique(total$Type.Name)))

# collist=c("#524FA1", "#FDB913", "lightgreen", "#00ADDC", "Darkgreen", "#ED1C24",
#           "Black", "Pink", "Brown", "Cyan", "midnightblue", "palevioletred3",
#           "lightcoral", "yellow4", "wheat4")


# save the tree data
write.tree(tr, file=tree_file)

# calcualte the plot size
numnodes <- length(m[,1])
plotsize = max(c(numnodes/20, 7))

# plot the tree
pdf(file=pdf_file, width = plotsize, height = plotsize)
plot(tr, type='fan', cex=0.8) #, tip.color=col)
dev.off()
