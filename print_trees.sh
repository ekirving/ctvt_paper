#!/usr/bin/env bash

for tree in trees/*.tre*
do
    echo "$tree"
    Rscript rscript/plot-phylo-tree-manual.r \
        njtree/all-pops.merged_v2.geno.data \
        phylogram \
        AndeanFox \
        "${tree}" \
        "${tree}.pdf"
done