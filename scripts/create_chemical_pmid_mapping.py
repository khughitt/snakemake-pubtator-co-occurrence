"""
Create a mapping from chemical mesh ids to pubmed article ids
"""
import json
import numpy as np
import pandas as pd

# load chemical data
chemical_dat = pd.read_feather(snakemake.input[0])

# convert concept id to pyarrow/string type (~7x faster in testing..)
chemical_dat.concept_id = chemical_dat.concept_id.astype('string[pyarrow]')

# baseline
# %timeit chemical_dat.concept_id == mesh_id
# 3.13 s ± 49.2 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

# string
# %timeit chemical_dat.concept_id == mesh_id
# 4.25 s ± 33.5 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

# string[pyarrow]
# %timeit  chemical_dat.concept_id == mesh_id
# 428 ms ± 48.3 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

# iterate over chemicals
mesh_ids = list(chemical_dat.concept_id.unique())
num_chemicals = len(mesh_ids)

# iterate over chemicals and retrieve associated pubmed ids for each
chemical_pmids = {}

for i, mesh_id in enumerate(mesh_ids):
    print(f"Processing MeSH ID {i + 1}/{len(mesh_ids)}")

    # slow step..
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
