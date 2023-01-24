"""
Filters Pubtator bioconcepts file
"""
import numpy as np
import pandas as pd

# re-assign "snakemake" variable to hide warnings
snek = snakemake

# load data
#  df = pd.read_csv(snek.input[0], sep='\t',
df = pd.read_csv("/data/raw/pubtator/2022-11-03/bioconcepts2pubtatorcentral", sep='\t',
                 names=['pmid', 'type', 'concept_id', 'mentions', 'resource'],
                 dtype={'pmid': np.int32, 'type': 'category', 'concept_id': 'category',
                        'mentions': str, 'resource': 'category'})

# if debug mode is enabled, sub-sample the data to speed up development
if snek.config['debug']:
    print(f"[DEBUG] Randomly sampling {snek.config['sample_pmids']} pubmed articles...")

    sampled_pmids = np.random.choice(df.pmid.unique(), snek.config['sample_pmids'])
    mask = df.pmid.isin(sampled_pmids)

    df = df[mask]

# keep track of original dataset size
num_pmids_orig = df.pmid.nunique()
num_entries_orig = df.shape[0]

# remove entries whose "mentions" field contains tabs/newline characters;
# these are often associated with erroneous entries
mask = ~(df.mentions.str.contains("\n") | df.mentions.str.contains("\t"))

num_dropped = (~mask).sum()
pct_dropped = 100 * (num_dropped / df.shape[0])

MSG = "Excluding %d / %d (%0.2f%%) entries whose 'mentions' field contains tabs/newline characters."
print(MSG % (num_dropped, df.shape[0], pct_dropped))

df = df[mask]

# limit analysis to articles relating to human research
# note: since it is possible an article could include annotations for both
# "human" and "mouse", etc., it may be useful to limit the analysis to articles
# _only_ including a reference to "human"..
human_pmids = set(df.pmid[df.concept_id == "9606"])
mask = df.pmid.isin(human_pmids)

num_dropped = num_pmids_orig - len(human_pmids)

pct_dropped = 100 * (num_dropped / num_pmids_orig)

print("Excluding %d / %d (%0.2f%%) pmids not relating to human research." %
        (num_dropped, num_pmids_orig, pct_dropped))

df = df[mask]

# filter out any pubmed articles with a greater than expected number of concept
# ids assigned to them
pmid_counts = df.pmid.value_counts()

to_keep = pmid_counts.index[pmid_counts <= snek.config['pmid_max_concepts']]

mask = df.pmid.isin(to_keep)

# article filtering stats
num_kept = len(to_keep)
num_total = df.pmid.nunique()
num_dropped = num_total - num_kept

pct_dropped = 100 * (num_dropped / num_total)

MSG = "Excluding %d / %d (%0.2f%%) pubmed articles with > %d concepts associated with them.."
print(MSG % (num_dropped, num_total, pct_dropped, snek.config['pmid_max_concepts']))

df = df.loc[mask]

# next, drop concepts which appear only a small number of times
concept_counts = df.concept_id.value_counts()

to_keep = concept_counts.index[concept_counts >= snek.config['concept_id_min_freq']]

mask = df.concept_id.isin(to_keep)

# concept_id filtering stats
mask_counts = mask.value_counts()

num_kept = 0 if True not in mask_counts else mask_counts[True]
num_total = df.shape[0]

num_dropped = 0 if False not in mask_counts else mask_counts[False]

pct_dropped = 100 * (num_dropped / num_total)

MSG = "Excluding %d / %d (%0.2f%%) concept ids which appear <%d times.." 
print(MSG % (num_dropped, num_total, pct_dropped, snek.config['concept_id_min_freq']))

df = df[mask]

# next, remove entries associated with excessively long "mentions" fields
mask = df.mentions.str.len() <= snek.config['mentions_max_length']

num_kept = mask.sum()
num_dropped = (~mask).sum()
num_total = df.shape[0]

pct_dropped = 100 * (num_dropped / num_total)

MSG = "Excluding %d / %d (%0.2f%%) entries with 'mentions' fields > %d characters.."
print(MSG % (num_dropped, num_total, pct_dropped, snek.config['mentions_max_length']))

df = df[mask]

# print summary of filtering results
print(f"Final dataset size: {df.shape[0]} rows ({len(set(df.pmid))} articles)")
print(f"Removed {num_entries_orig} rows in total ({num_pmids_orig} articles)") 

# save filtered dataset
#df.reset_index(drop=True).to_csv(snek.output[0], sep='\t')
#df.to_csv(snek.output[0], index=False, sep='\t')
df.to_feather(snek.output[0])
