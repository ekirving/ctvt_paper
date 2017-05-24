#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
data_file=args[1]
pdf_file=args[2]
max_k=args[3]

# TODO remove when done testing
# setwd("/Users/Evan/Dropbox/Code/socotra")

# data_file = "admix/all-pops.BovineHD.pruned.sorted.3.data"
# pdf_file = "pdf/all-pops.BovineHD.admix.K.3.pdf"

# data_file = "admix/modern-pops.BovineHD.pruned.sorted.22.data"
# pdf_file = "pdf/modern-pops.BovineHD.admix.K.22.pdf"

# max_k = 22

# read the data file
dat = read.table(data_file, header=TRUE)

# sum the components for each population
sums <- aggregate(. ~ Population, dat[-2], mean)
rownames(sums) <- sums[,1]

# get the largest k component for each population
maxk <- as.data.frame(apply(sums[-1],1, which.max))
colnames(maxk) <- c('max.k')

maxval <- as.data.frame(apply(sums[-1],1, max))
colnames(maxval) <- c('max.val')

# merge those columns
dat <- merge(dat, maxk, by.x = 'Population', by.y = 'row.names')
dat <- merge(dat, maxval, by.x = 'Population', by.y = 'row.names')

# use the new column to sort the matrix
dat <- dat[order(dat$max.k, -rank(dat$max.val)),]

# drop the columns now we're done
dat <- dat[, !(colnames(dat) %in% c('max.k', 'max.val'))]

# # get the population metadata
# info <- read.table('list_pop_gps3.txt', sep = "\t", quote = '"', header=TRUE)
#
# # merge the metadata
# dat <- merge(dat, info[c(1,4)], by.x = 'Population', by.y = 'Code')
#
# # sort by the population groups
# dat <- dat[order(dat$Type.Name),]
#
# # drop the exta column we added
# dat <- dat[, !(colnames(dat) %in% c('Type.Name'))]

# perserve the existing ordering
dat$Sample <- factor(dat$Sample, levels = dat$Sample)
dat$Population <- factor(dat$Population, levels = unique(dat$Population))

dat.m <- melt(dat,id=c("Population", "Sample"))

if (strtoi(max_k) <= 12) {
    # use the nice 12-colour palette from colorbrewer2
    # http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12
    colours=c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00',
              '#cab2d6', '#6a3d9a', '#ffff99', '#b15928')

} else {
    # use a large 104-colour palette from Ssnot! (w/ minor modifications)
    # http://godsnotwheregodsnot.blogspot.co.uk/2012/09/color-distribution-methodology.html
    colours=c(           "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                         "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
              "#5A0007", "#809693",                       "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
              "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
              "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
              "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
              "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
              "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
              "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
              "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
              "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
              "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
              "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
              '#000000', "#FFDBE5")
}

# how many rows do we need
numpops <- length(levels(dat[,1]))
nrows <- ceiling(numpops/20)

pdf(file=pdf_file, width = 16, height = 9)

ggplot(dat.m, aes(x = Sample, y = value, fill=variable, group=Population)) +

    # fixed width, merged bars
    # geom_bar(stat='identity', width=1) +
    # facet_grid(~Population, scales = "free") +
    # facet_grid(~Population, scales = "free") +

    # variable width, merged bars
    # geom_bar(stat='identity', width=1) +
    # facet_grid(~Population, scales = "free", space = "free") +
    # facet_grid(~Population, scales = "free", space = "free") +

    # variable width, seperate bars
    # geom_bar(stat='identity') +
    # facet_grid(~Population, scales = "free", space = "free") +

    # fixed width, merged bars, wrapped across n rows
    geom_bar(stat='identity', width=1) +
    facet_wrap(~Population, scales = "free", nrow = nrows) +

    # fixed width, seperate bars, wrapped across n rows
    # geom_bar(stat='identity') +
    # facet_wrap(~Population, scales = "free", nrow = nrows) +

    theme(
          # legend.title=element_blank(),
          # legend.key = element_blank(),
          legend.position = "none",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          panel.spacing = unit(0, "lines")
          ) +

    xlab(paste("K =", dim(dat)[2]-2)) +
    ylab("") +
    # ylab("Ancestry") +

    scale_fill_manual(values=colours)

dev.off()
