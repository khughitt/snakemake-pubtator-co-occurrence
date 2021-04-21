"""
Pubtator Co-occurrence matrix pipeline (v2)

KH (April 3, 2021)

Constructs a human gene-gene co-occurence matrix from pubtator central bioconcepts data.
"""
import itertools
import os
import random
import time
import numpy as np
import pandas as pd
from os.path import join
from scipy import sparse
from scipy.io import mmread, mmwrite
from sklearn.preprocessing import scale

# os.environ["MODIN_ENGINE"] = "dask"
# import modin.pandas as pd

random.seed(config['random_seed'])

out_dir = join(config["out_dir"], config["version"])

#
# Snakemake rules
#
rule all:
    input:
        join(out_dir, 'filtered', 'bioconcepts2pubtatorcentral_filtered.feather')
        # expand(join(out_dir, "submats", "{topic1}_{topic2}.mtx"), topic1=topics, topic2=topics)
        # join(out_dir, "submats", "{topic1}_{topic2}.mtx")

# rule split_dataset:
#     input:
#         join(out_dir, 'bioconcepts2pubtatorcentral_filtered.feather'),
#         join(out_dir, 'bioconcepts2pubtatorcentral_article_clusters.feather')
#     output: expand(join(out_dir, 'subsets', 'bioconcepts2pubtatorcentral_filtered_{i}.feather'), i=range(config['clustering']['num_clusters']))
#     run:
#         # load full dataset
#         dat = pd.read_feather(input[0]).drop(['index'], axis=1)
#
#         pmids = list(dat.pmid.unique())
#
#         # load paper cluster assignments
#         clusters = pd.read_feather(input[1]).set_index('pmid')
#
#         # split into n chunks with similar numbers of articles in each
#         for i in clusters.cluster.unique():
#             # get articles in cluster and store result
#             subset_ids = clusters[clusters.cluster == i].index
#
#             # cluster indices start from 0 / snakemake starts from 1
#             dat[dat.pmid.isin(subset_ids)].reset_index().to_feather(output[i])

rule filter_dataset:
    input:
        config['pubtator_data_path']
    output:
        join(out_dir, 'filtered', 'bioconcepts2pubtatorcentral_filtered.feather')
        # join(out_dir, 'filtered', 'bioconcepts2pubtatorcentral_vocab.feather')
    run:
        # load data;
        # note: a very small number of malformed-lines containing >4 tab separators
        # were found when processing the data and are explicitly ignored below
        dat = pd.read_csv(input[0], sep='\t', error_bad_lines=False,
                          names=['pmid', 'type', 'concept_id', 'mentions', 'resource'],
                          dtype={'pmid': np.int32, 'type': str, 'concept_id': str,
                                 'mentions': str, 'resource': str})

        # limit analysis to articles relating to human research
        # note: since it is possible an article could include annotations for both
        # "human" and "mouse", etc., it may be useful to limit the analysis to articles
        # _only_ including a reference to "human"..
        human_pmids = set(dat.pmid[dat.concept_id == "9606"])
        mask = dat.pmid.isin(human_pmids)

        num_total = dat.pmid.nunique()
        num_dropped = num_total - len(human_pmids)

        pct_dropped = 100 * (num_dropped / num_total)

        print("Excluding %d / %d (%0.2f%%) pmids not relating to human research." % (num_dropped, num_total, pct_dropped))

        dat = dat[mask]

        # filter out any pubmed articles with a greater than expected number of concept
        # ids assigned to them
        pmid_counts = dat.pmid.value_counts()

        to_keep = pmid_counts.index[pmid_counts <= config['pmid_max_concepts']]

        mask = dat.pmid.isin(to_keep)

        # article filtering stats
        num_kept = len(to_keep)
        num_total = dat.pmid.nunique()
        num_dropped = num_total - num_kept

        pct_dropped = 100 * (num_dropped / num_total)

        print("Excluding %d / %d (%0.2f%%) pubmed articles with > %d concepts associated with them.." % (num_dropped, num_total, pct_dropped, config['pmid_max_concepts']))

        dat = dat.loc[mask]

        # next, drop concepts which appear only a small number of times
        concept_counts = dat.concept_id.value_counts()

        to_keep = concept_counts.index[concept_counts >= config['concept_id_min_freq']]

        mask = dat.concept_id.isin(to_keep)

        # concept_id filtering stats
        mask_counts = mask.value_counts()

        num_kept = mask_counts[True]
        num_total = dat.shape[0]
        num_dropped = mask_counts[False]

        pct_dropped = 100 * (num_dropped / num_total)

        print("Excluding %d / %d (%0.2f%%) concept ids which appear <%d times.." % (num_dropped, num_total, pct_dropped, config['concept_id_min_freq']))

        dat = dat[mask]

        # TODO: print removal stats
        print(f"Final dataset size: {dat.shape[0]} rows ({len(set(dat.pmid))} articles)")

        # get a list of all unique concept ids remaining
        # vocab = sorted(dat.concept_id.unique())

        # map from concept ids to numeric indices to use in the output matrix
        # mapping = {}
        #
        # for i, annot in enumerate(vocab):
        #     mapping[annot] = i
        #
        # dat['ind'] = dat.concept_id.map(mapping)

        # sort indices within each article
        #
        # from the docs for sparse.lil_matrix:
        #
        # "This is a structure for constructing sparse matrices incrementally.
        # Note that inserting a single item can take linear time in the worst
        # case; to construct a matrix efficiently, make sure the items are
        # pre-sorted by index, per row."
        #
        # dat = dat.sort_values(['pmid', 'ind'])

        # save filtered dataset
        dat.reset_index(drop=True).to_feather(output[0])

        # also save a complete list of the vocabulary
        # vocab_df = dat[['type', 'concept_id', 'ind']].drop_duplicates().sort_values('ind')

        # vocab_df.reset_index().drop('index', axis=1).to_feather(output[1])

