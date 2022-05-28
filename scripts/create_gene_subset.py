"""
Creates a filtered version of the bioconcepts file containing only gene entries
"""
import pandas as pd

# load filtered dataset
dat = pd.read_csv(snakemake.input[0], sep='\t')

# load list of human entrez gene identifiers
with open(snakemake.input[1]) as fp:
    entrez_ids = [x.strip() for x in fp.readlines()]

# create a subset with only gene entries, restricted to GNormPlus/Human Entrez
# gene identifiers..
entrez_mask = dat.concept_id.isin(entrez_ids)

dat = dat[(dat.type == 'Gene') & (dat.resource == 'GNormPlus') & entrez_mask]

# store gene subset
dat.reset_index(drop=True).to_feather(snakemake.output[0])
