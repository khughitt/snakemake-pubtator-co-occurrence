"""
Pubtator Co-occurrence matrix pipeline

KH 2020.09.01

Processes abstracts and annotations retrieved from PubTator and generates a sparse
co-occurrence matrix representation of the article annotations contained in them.
"""
import itertools
import random
import time
import numpy as np
import pandas as pd
from os.path import join
from scipy import sparse
from sklearn.cluster import MiniBatchKMeans

random.seed(config['random_seed'])

rule combine_matrices:
    input:
        expand(join(config["out_dir"], 'co-occurrence', 'bioconcepts2pubtatorcentral_comat_{i}.npz'),
               i=range(config['clustering']['num_clusters'])),
        join(config["out_dir"], 'bioconcepts2pubtatorcentral_vocab.feather')

rule filter_dataset:
    input:
        config['pubtator_data_path']
    output:
        join(config["out_dir"], 'bioconcepts2pubtatorcentral_filtered.feather'),
        join(config["out_dir"], 'bioconcepts2pubtatorcentral_vocab.feather')
    run:
        # load data;
        # note: a very small number of malformed-lines containing >4 tab separators 
        # were found when processing the data and are explicitly ignored below
        dat = pd.read_csv(input[0], sep='\t', error_bad_lines=False,
                          names=['pmid', 'type', 'concept_id', 'mentions', 'resource'],
                          dtype={'pmid': np.int32, 'type': str, 'concept_id': str,
                                 'mentions': str, 'resource': str})

        # first, count the number of occurrences of each concept id
        concept_counts = dat.concept_id.value_counts()

        # drop concepts that only appear once; this will result in a ~80% reduction in
        # the number of terms to compare and result in a significant reduction in time.
        to_keep = concept_counts.index[concept_counts > 1]
        dat = dat[dat.concept_id.isin(to_keep)]

        # for now, exclude chemical / species concepts to reduce the problem size further
        dat = dat[np.logical_not(dat.type.isin(config['exclude_concepts']))]

        # get a list of all unique concept ids remaining
        vocab = sorted(dat.concept_id.unique())

        # map from concept ids to numeric indices to use in the output matrix
        mapping = {}

        for i, annot in enumerate(vocab):
            mapping[annot] = i

        dat['ind'] = dat.concept_id.map(mapping)

        # sort indices within each article
        #
        # from the docs for sparse.lil_matrix:
        #
        # "This is a structure for constructing sparse matrices incrementally.
        # Note that inserting a single item can take linear time in the worst
        # case; to construct a matrix efficiently, make sure the items are
        # pre-sorted by index, per row."
        #
        dat = dat.sort_values(['pmid', 'ind'])

        # save filtered dataset
        dat.reset_index().to_feather(output[0])

        # also save a complete list of the vocabulary
        vocab_df = dat[['concept_id', 'ind']].drop_duplicates().sort_values('ind')

        vocab_df.reset_index().to_feather(output[1])

rule cluster_articles:
    input:
        join(config["out_dir"], 'bioconcepts2pubtatorcentral_filtered.feather')
    output:
        join(config["out_dir"], 'bioconcepts2pubtatorcentral_article_matrix.feather'),
        join(config["out_dir"], 'bioconcepts2pubtatorcentral_article_clusters.feather')
    run:
        # load data
        dat = pd.read_feather(input[0]).drop(['index'], axis=1)

        #
        # rank concepts by frequency
        #
        # for the top n = 1000 concepts, this produces a subset with ~56% (81259480)
        # rows; in other words, a little over half of the annotations come from the same
        # 1000 concepts.
        #
        concept_counts = dat.ind.value_counts().sort_values(ascending = False)


        # TODO: how about using top N PC's instead of the concepts themselves to
        # cluster?
        
        # subset data to include only the top N concepts
        top_concepts = concept_counts.head(config['clustering']['num_concepts']).index
        dat_subset = dat[dat.ind.isin(top_concepts)]

        # generate a binary article x concept matrix
        dat_subset['present'] = 1

        # NOTE: there is a limitation in pandas which prevents pivoting from working
        # for large dataframes like this one.
        # 
        # see: https://github.com/pandas-dev/pandas/issues/26314
        X = dat_subset.pivot(index='pmid', columns='ind', values=['present'])

        # store pubmed article order
        pmids = X.index

        # set colnames to string versions of concept ids
        X.columns = [str(x[1]) for x in X.columns]

        # store article x concept dataframe
        X.reset_index().to_feather(output[0])

        # replace nan's with zeros
        X[X.isna()] = 0

        # convert to numpy array
        X = X.to_numpy()

        # NOTE: for the ~13,000,000 x 160 matrix, the sparsity is ~98%.
        kmeans = MiniBatchKMeans(n_clusters=config['clustering']['num_clusters'], random_state=0, batch_size=1000000, max_iter=10).fit(X)

        # store cluster labels
        clusters = pd.DataFrame({'cluster': kmeans.labels_}, index=pmids)

        clusters.reset_index().to_feather(output[1])

rule split_dataset:
    input:
        join(config["out_dir"], 'bioconcepts2pubtatorcentral_filtered.feather'),
        join(config["out_dir"], 'bioconcepts2pubtatorcentral_article_clusters.feather')
    output: expand(join(config["out_dir"], 'bioconcepts2pubtatorcentral_filtered_{i}.feather'), i=range(config['clustering']['num_clusters'])) 
    run:
        # load full dataset
        dat = pd.read_feather(input[0]).drop(['index'], axis=1)

        pmids = list(dat.pmid.unique())

        # load paper cluster assignments
        clusters = pd.read_feather(input[1]).set_index('pmid')

        # for i, subset_ids in enumerate(pmid_chunks):

        # split into n chunks with similar numbers of articles in each
        for i in clusters.cluster.unique():
            # get articles in cluster and store result
            subset_ids = clusters[clusters.cluster == i].index

            # cluster indices start from 0 / snakemake starts from 1
            dat[dat.pmid.isin(subset_ids)].reset_index().to_feather(output[i])

rule build_cooccurrence_matrices:
    input:
        join(config["out_dir"], 'bioconcepts2pubtatorcentral_filtered_{i}.feather'),
        join(config["out_dir"], 'bioconcepts2pubtatorcentral_vocab.feather')
    output:
        join(config["out_dir"], 'co-occurrence', 'bioconcepts2pubtatorcentral_comat_{i}.npz') 
    run:
        # load data subset
        dat = pd.read_feather(input[0]).drop(['index'], axis=1)

        # load complete vocabulary
        vocab_df = pd.read_feather(input[1]).drop(['index'], axis=1)

        vocab = vocab_df.concept_id
        vocab_size = len(vocab)

        # create a lil sparse matrix
        mat = sparse.lil_matrix((vocab_size, vocab_size), dtype=np.int32)

        # get a list of unique pubmed ids
        pmids = dat.pmid.unique()
        num_articles = len(pmids)

        # iterate over articles and populate matrix
        t0 = time.time()
        
        for counter, pmid in enumerate(pmids):
            if counter % 10000 == 0:
                print("Processing article %d / %d (%0.2f%%)..." % (counter + 1,
                    num_articles, 100 * (counter + 1) / num_articles))

            annots = dat[dat.pmid == pmid].ind.unique()

            # for i, j in itertools.combinations(range(len(annots)), 2):
            for i, j in itertools.combinations(annots, 2):
                mat[i, j] += 1 

        t1 = time.time()
        tdelt = (t1 - t0) * 60 * 60
        rate = num_articles / (t1 - t0)

        print("Finished building co-occurrence matrix.. (Time: %0.1f hours, Rate: %0.2f articles/second)" % (tdelt, rate))

        # convert to coo format and store result
        sparse.save_npz(output[0], mat.tocoo())

