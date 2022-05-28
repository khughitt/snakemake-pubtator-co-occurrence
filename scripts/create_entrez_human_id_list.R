#!/bin/env Rscript
#
# Generates a list of Human Entrez gene identifers
#
library(org.Hs.eg.db)
entrez = keys(org.Hs.eg.db, keytype='ENTREZID')
writeLines(entrez, snakemake@output[[1]])
