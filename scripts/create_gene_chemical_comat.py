"""
Creates a gene-chemical co-occurrence matrix
"""
import json
import pandas as pd
import numpy as np

# load gene and drug pmid mappings
with open(snakemake.input[0]) as fp:
    all_gene_pmids = json.load(fp)

with open(snakemake.input[1]) as fp:
    all_chemical_pmids = json.load(fp)

# create empty matrix to store gene-chemical co-occurrence counts
entrez_ids = all_gene_pmids.keys()
num_genes = len(entrez_ids)

mesh_ids = all_chemical_pmids.keys()

# exclude "-" and chemicals with fewer associated pmids
# TODO: move upstream prior to next major run..

# note: in current version of pipeline, setting a minimum pmid count of "50" for
# the chemical data results in 17,811/77,041 chemicals being included.. 
chem_pmid_counts = pd.Series([len(a) for a in all_chemical_pmids.values()],
                                index=all_chemical_pmids.keys())

MIN_PMID_COUNT = snakemake.config['chem_comat_min_pmid']

to_keep = chem_pmid_counts[chem_pmid_counts >= MIN_PMID_COUNT].index
mesh_ids = [x for x in mesh_ids if x in to_keep]

mesh_ids.remove('-')

num_chemicals = len(mesh_ids)

comat = np.empty((num_genes, num_chemicals))
comat.fill(np.nan)

print(f"Processing {num_genes} genes...")

# iterate over pairs of genes
for i, gene in enumerate(entrez_ids):
    # get pubmed ids associated with gene
    gene_pmids = all_gene_pmids[gene]

    if i % 100 == 0:
        print(f"gene {i}/{num_genes}...")

    for j, chemical in enumerate(mesh_ids):
        # get pubmed ids associated with chemical
        chemical_pmids = all_chemical_pmids[chemical]

        # compute gene-chemical co-occurrence count
        num_shared = len(set(gene_pmids).intersection(chemical_pmids))

        comat[i, j] = num_shared

# store gene-chemical co-occurrence matrix
comat = pd.DataFrame(comat, index=entrez_ids, columns=mesh_ids)
comat.reset_index().rename(columns={'index': 'entrez_id'}).to_feather(output[0])
