#!/usr/bin/env Rscript
library(combinat)
library(foreach)
library(doParallel)

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
qfile=args[1]
kprior=args[2]
qsorted=args[3]

# TODO remove when done testing
# setwd("/Users/Evan/Dropbox/Code/ctvt")
# qfile = "admix/qpgraph-pops.merged_v2_hq2_nomex_ctvt.pruned.4.Q"
# kprior = "admix/qpgraph-pops.merged_v2_hq2_nomex_ctvt.pruned.sorted.3.Q"
# qsorted = "admix/qpgraph-pops.merged_v2_hq2_nomex_ctvt.pruned.sorted.4.Q"

# read the Q matrices into dataframes
dat1 = read.table(kprior)
dat2 = read.table(qfile)

# pad the smaller matrix with zeros (to allow simple matrix algebra)
numPad = length(dat2) - length(dat1)
for (i in 1:numPad) {
  dat1[sprintf("P%d",i)] <- 0
}

# if K < 10 we can test all the permutations to find the absolute best fit
# this algorithm has complexity of O(n!) (e.g. 10! = 3,628,800 perms)
if (length(dat2) < 10) {

  # make a vector containing all the possible permutations of the column orders
  colPerms <- permn(names(dat2))

  sums <- c()

  # initialise a parallel cluster
  cl <- makeCluster(detectCores()[1] * 0.5) # only use 50% of the available cores
  registerDoParallel(cl)

  # TODO this would be a lot quicker if the job was chunked
  # iterate over every permutation of the matrix in parallel
  sums <- foreach(i = seq_along(colPerms), .combine=c) %dopar% {
    # calculate the sum of the absolute difference
    sum(abs(dat1 - dat2[colPerms[[i]]]))
  }

  stopCluster(cl)

  # find the index with the lowest overall difference
  idx <- which(sums == min(sums), arr.ind = TRUE)[1]

  # apply the optimal column order
  sorted2 <- dat2[colPerms[[idx]]]

} else {

  # using the k-1 matrix, find the best fit for each column in the k matrix (without replacement)
  # this simple algorithm optimises the order on a per column basis, rather than for the matrix as a whole, i.e. O(n)
  tmp <- dat2
  best <- c()

  for (i in names(dat1)) {
    # subtract the current column from the tmp matrix
    sums <- colSums(abs(tmp-dat1[[i]]))

    # find the best column match
    idx <- which(sums == min(sums), arr.ind = FALSE)
    best[[i]] <- names(idx)

    # delete the column from the matrix
    tmp <- tmp[-idx]
  }

  # apply the optimal column order
  sorted2 <- dat2[best]
}

# write the matrix back to the filesystem, in the optimal order
write.table(sorted2, qsorted, row.names = FALSE, col.names = FALSE)
