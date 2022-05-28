"""
Counts the number of unique genes associated with each pmid, after filtering
"""
import pandas as pd

dat = pd.read_feather(snakemake.input[0])

counts = dat.groupby('pmid').concept_id.nunique()
counts.rename("num").reset_index().to_feather(snakemake.output[0])
