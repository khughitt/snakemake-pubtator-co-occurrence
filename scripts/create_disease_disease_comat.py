"""
Creates a disease-disease co-occurrence matrix
"""
import json
import pandas as pd
import numpy as np

# load disease/pmid mapping
with open(snakemake.input[0]) as fp:
    disease_pmids = json.load(fp)

# create empty matrix to store disease-disease co-occurrence counts
mesh_ids = disease_pmids.keys()
num_diseases = len(mesh_ids)

comat = np.empty((num_diseases, num_diseases))
comat.fill(np.nan)

# iterate over pairs of diseases
for i, disease1 in enumerate(mesh_ids):
    # get pubmed ids associated with disease 1
    disease1_pmids = disease_pmids[disease1]

    for j, disease2 in enumerate(mesh_ids):
        # skip symmetric comparisons
        if not np.isnan(comat[i, j]):
            continue

        disease2_pmids = disease_pmids[disease2]

        # compute disease-disease co-occurrence count
        num_shared = len(set(disease1_pmids).intersection(disease2_pmids))

        comat[i, j] = comat[j, i] = num_shared

# store disease-disease co-occurrence matrix
comat = pd.DataFrame(comat, index=mesh_ids, columns=mesh_ids)
comat.reset_index().rename(columns={'index': 'mesh_id'}).to_feather(snakemake.output[0])
