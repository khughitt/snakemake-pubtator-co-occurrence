#!/bin/env R
#
# Generates a mapping from disease MeSH ID's to the most frequent associated mentions
# associated with those ids.
#
# In the case of ties, an arbitrary mention is selected.
#
# KH (June 2021)
#
library(arrow)
library(tidyverse)

dat <- read_feather(snakemake@input[[1]]) %>%
  select(concept_id, mentions)

# exclude entries with missing "mentions"
dat <- dat[!is.na(dat$mentions), ]

# exclude entries with no assigned mesh id
mask <- grepl("MESH:", dat$concept_id)
dat <- dat[mask, ]

# choose the most common mention for each mesh term
mapping <- dat %>%
  group_by(concept_id, mentions) %>%
  summarize(n=n()) %>%
  arrange(desc(n)) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  select(concept_id, mentions)

colnames(mapping) <- c("mesh_id", snakemake@params$field)

write_tsv(mapping, snakemake@output[[1]])
