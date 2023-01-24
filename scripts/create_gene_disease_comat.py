"""
Creates a gene-disease co-occurrence matrix
"""
import json
import pandas as pd
import numpy as np

# load gene and drug pmid mappings
with open(snakemake.input[0]) as fp:
    all_gene_pmids = json.load(fp)

with open(snakemake.input[1]) as fp:
    all_disease_pmids = json.load(fp)

# create empty matrix to store gene-disease co-occurrence counts
entrez_ids = all_gene_pmids.keys()
num_genes = len(entrez_ids)

mesh_ids = all_disease_pmids.keys()
num_diseases = len(mesh_ids)

comat = np.empty((num_genes, num_diseases), dtype=np.uint32)
comat.fill(np.nan)

# iterate over pairs of genes
for i, gene in enumerate(entrez_ids):
    # get pubmed ids associated with gene
    gene_pmids = all_gene_pmids[gene]

    for j, disease in enumerate(mesh_ids):
        # get pubmed ids associated with disease
        disease_pmids = all_disease_pmids[disease]

        # compute gene-disease co-occurrence count
        num_shared = len(set(gene_pmids).intersection(disease_pmids))

        comat[i, j] = num_shared

# store gene-disease co-occurrence matrix
comat = pd.DataFrame(comat, index=entrez_ids, columns=mesh_ids)
comat.reset_index().rename(columns={'index': 'entrez_id'}).to_feather(snakemake.output[0])
