"""
Creates a chemical-chemical co-occurrence matrix
"""
import json
import pandas as pd
import numpy as np

# load chemical/pmid mapping
with open(snakemake.input[0]) as fp:
    unfiltered_chemical_pmids = json.load(fp)

# filter any chemicals with less than N citations; helps keep the co-occurrence
# matrix size more manageable while only losing relatively low-information
# low-citation counts..
chemical_pmids = {}

for mesh_id in unfiltered_chemical_pmids:
    if len(unfiltered_chemical_pmids[mesh_id]) > snakemake.config['concept_id_min_freq_chem']:
        chemical_pmids[mesh_id] = unfiltered_chemical_pmids[mesh_id]

# create empty matrix to store chemical-chemical co-occurrence counts
mesh_ids = chemical_pmids.keys()
num_chemicals = len(mesh_ids)

comat = np.empty((num_chemicals, num_chemicals))
comat.fill(np.nan)

# iterate over pairs of chemicals
for i, chemical1 in enumerate(mesh_ids):
    # get pubmed ids associated with chemical 1
    chemical1_pmids = chemical_pmids[chemical1]

    for j, chemical2 in enumerate(mesh_ids):
        # skip symmetric comparisons
        if not np.isnan(comat[i, j]):
            continue

        chemical2_pmids = chemical_pmids[chemical2]

        # compute chemical-chemical co-occurrence count
        num_shared = len(set(chemical1_pmids).intersection(chemical2_pmids))

        comat[i, j] = comat[j, i] = num_shared

# store chemical-chemical co-occurrence matrix
comat = pd.DataFrame(comat, index=mesh_ids, columns=mesh_ids)
comat.reset_index().rename(columns={'index': 'mesh_id'}).to_feather(snakemake.output[0])
