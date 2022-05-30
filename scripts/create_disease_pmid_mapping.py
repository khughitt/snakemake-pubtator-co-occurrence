"""
Create a mapping from disease mesh ids to pubmed article ids
"""
import json
import numpy as np
import pandas as pd

# load disease data
disease_dat = pd.read_feather(snakemake.input[0])

# convert concept id to pyarrow/string type (~7x faster in testing..)
disease_dat.concept_id = disease_dat.concept_id.astype('string[pyarrow]')

# iterate over diseases
mesh_ids = list(disease_dat.concept_id.unique())
num_diseases = len(mesh_ids)

# iterate over diseases and retrieve associated pubmed ids for each
disease_pmids = {}

for mesh_id in mesh_ids:
    mask = disease_dat.concept_id == mesh_id
    disease_pmids[mesh_id] = set(disease_dat[mask].pmid.values)

# encoder to convert int64 elements to generic ints and sets to lists during
# json serialization
def encoder(object):
    if isinstance(object, np.generic):
        return object.item()
    elif isinstance(object, set):
        return list(object)

# store disease -> pmid mapping as json
with open(snakemake.output[0], "w") as fp:
    fp.write(json.dumps(disease_pmids, default=encoder))
