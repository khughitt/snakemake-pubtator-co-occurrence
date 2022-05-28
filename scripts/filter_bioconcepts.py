"""
Filters Pubtator bioconcepts file
"""
import numpy as np
import pandas as pd

# load data
dat = pd.read_csv(snakemake.input[0], sep='\t',
                  names=['pmid', 'type', 'concept_id', 'mentions', 'resource'],
                  dtype={'pmid': np.int32, 'type': str, 'concept_id': str,
                         'mentions': str, 'resource': str})

# if debug mode is enabled, sub-sample the data to speed up development
if snakemake.config['debug']:
    print(f"[DEBUG] Randomly sampling {snakemake.config['sample_pmids']} pubmed articles...")

    sampled_pmids = np.random.choice(dat.pmid.unique(), snakemake.config['sample_pmids'])
    mask = dat.pmid.isin(sampled_pmids)

    dat = dat[mask]

# keep track of original dataset size
num_pmids_orig = dat.pmid.nunique()
num_entries_orig = dat.shape[0]

# limit analysis to articles relating to human research
# note: since it is possible an article could include annotations for both
# "human" and "mouse", etc., it may be useful to limit the analysis to articles
# _only_ including a reference to "human"..
human_pmids = set(dat.pmid[dat.concept_id == "9606"])
mask = dat.pmid.isin(human_pmids)

num_dropped = num_pmids_orig - len(human_pmids)

pct_dropped = 100 * (num_dropped / num_pmids_orig)

print("Excluding %d / %d (%0.2f%%) pmids not relating to human research." %
        (num_dropped, num_pmids_orig, pct_dropped))

dat = dat[mask]

# filter out any pubmed articles with a greater than expected number of concept
# ids assigned to them
pmid_counts = dat.pmid.value_counts()

to_keep = pmid_counts.index[pmid_counts <= snakemake.config['pmid_max_concepts']]

mask = dat.pmid.isin(to_keep)

# article filtering stats
num_kept = len(to_keep)
num_total = dat.pmid.nunique()
num_dropped = num_total - num_kept

pct_dropped = 100 * (num_dropped / num_total)

msg = "Excluding %d / %d (%0.2f%%) pubmed articles with > %d concepts associated with them.."
print(msg % (num_dropped, num_total, pct_dropped, snakemake.config['pmid_max_concepts']))

dat = dat.loc[mask]

# next, drop concepts which appear only a small number of times
concept_counts = dat.concept_id.value_counts()

to_keep = concept_counts.index[concept_counts >= snakemake.config['concept_id_min_freq']]

mask = dat.concept_id.isin(to_keep)

# concept_id filtering stats
mask_counts = mask.value_counts()

num_kept = 0 if True not in mask_counts else mask_counts[True]
num_total = dat.shape[0]

num_dropped = 0 if False not in mask_counts else mask_counts[False]

pct_dropped = 100 * (num_dropped / num_total)

msg = "Excluding %d / %d (%0.2f%%) concept ids which appear <%d times.." 
print(msg % (num_dropped, num_total, pct_dropped, snakemake.config['concept_id_min_freq']))

dat = dat[mask]

# next, remove entries associated with excessively long "mentions" fields
mask = dat.mentions.str.len() <= snakemake.config['mentions_max_length']

num_kept = mask.sum()
num_dropped = (~mask).sum()
num_total = dat.shape[0]

pct_dropped = 100 * (num_dropped / num_total)

msg = "Excluding %d / %d (%0.2f%%) entries with 'mentions' fields > %d characters.."
print(msg % (num_dropped, num_total, pct_dropped, snakemake.config['mentions_max_length']))

dat = dat[mask]

# print summary of filtering results
print(f"Final dataset size: {dat.shape[0]} rows ({len(set(dat.pmid))} articles)")
print(f"Removed {num_entries_orig} rows in total ({num_pmids_orig} articles)") 

# save filtered dataset
#dat.reset_index(drop=True).to_csv(snakemake.output[0], sep='\t')
dat.to_csv(snakemake.output[0], index=False, sep='\t')
