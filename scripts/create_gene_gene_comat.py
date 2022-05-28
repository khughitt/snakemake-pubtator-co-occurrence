"""
Creates a gene-gene co-occurrence matrix
"""
import json
import pandas as pd
import numpy as np

# load gene/pmid mapping
with open(snakemake.input[0]) as fp:
    gene_pmids = json.load(fp)

# create empty matrix to store gene-gene co-occurrence counts
entrez_ids = gene_pmids.keys()
num_genes = len(entrez_ids)

comat = np.empty((num_genes, num_genes))
comat.fill(np.nan)

# iterate over pairs of genes
for i, gene1 in enumerate(entrez_ids):
    # get pubmed ids associated with gene 1
    gene1_pmids = gene_pmids[gene1]

    for j, gene2 in enumerate(entrez_ids):
        # skip symmetric comparisons
        if not np.isnan(comat[i, j]):
            continue

        gene2_pmids = gene_pmids[gene2]

        # compute gene-gene co-occurrence count
        num_shared = len(set(gene1_pmids).intersection(gene2_pmids))

        comat[i, j] = comat[j, i] = num_shared

# store gene-gene co-occurrence matrix
comat = pd.DataFrame(comat, index=entrez_ids, columns=entrez_ids)
comat.reset_index().rename(columns={'index': 'entrez_id'}).to_feather(snakemake.output[0])
