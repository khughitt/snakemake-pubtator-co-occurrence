"""
Create a mapping from entrez gene ids to pubmed article ids
"""
import json
import numpy as np
import pandas as pd

# load gene data
gene_dat = pd.read_feather(snakemake.input[0])

# iterate over genes
entrez_ids = list(gene_dat.concept_id.unique())
num_genes = len(entrez_ids)

# iterate over genes and retrieve associated pubmed ids for each
gene_pmids = {}

for entrez_id in entrez_ids:
    mask = gene_dat.concept_id == entrez_id
    gene_pmids[entrez_id] = set(gene_dat[mask].pmid.values)

# encoder to convert int64 elements to generic ints and sets to lists during
# json serialization
def encoder(object):
    if isinstance(object, np.generic):
        return object.item()
    elif isinstance(object, set):
        return list(object)

# store gene -> pmid mapping as json
with open(snakemake.output[0], "w") as fp:
    fp.write(json.dumps(gene_pmids, default=encoder))
