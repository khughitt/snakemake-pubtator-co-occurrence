"""
Create a mapping from chemical mesh ids to pubmed article ids
"""
import json
import numpy as np
import pandas as pd

# load chemical data
chemical_dat = pd.read_feather(snakemake.input[0])

# iterate over chemicals
mesh_ids = list(chemical_dat.concept_id.unique())
num_chemicals = len(mesh_ids)

# iterate over chemicals and retrieve associated pubmed ids for each
chemical_pmids = {}

for mesh_id in mesh_ids:
    mask = chemical_dat.concept_id == mesh_id
    chemical_pmids[mesh_id] = set(chemical_dat[mask].pmid.values)

# encoder to convert int64 elements to generic ints and sets to lists during
# json serialization
def encoder(object):
    if isinstance(object, np.generic):
        return object.item()
    elif isinstance(object, set):
        return list(object)

# store chemical -> pmid mapping as json
with open(snakemake.output[0], "w") as fp:
    fp.write(json.dumps(chemical_pmids, default=encoder))
