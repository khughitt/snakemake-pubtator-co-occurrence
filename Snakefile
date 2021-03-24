"""
Pubtator Co-occurrence matrix pipeline

KH 2020.09.07

Processes abstracts and annotations retrieved from PubTator and generates a sparse
co-occurrence matrix representation of the article annotations contained in them.
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
from sklearn.cluster import MiniBatchKMeans

random.seed(config['random_seed'])

out_dir = join(config["out_dir"], config["version"])

# get pairs of annotation topics (types)
topics = config['include_concepts'] 
topic_pairs = list(itertools.combinations(topics, 2))

# add self-comparisons (e.g. gene x gene)
for topic in topics:
    topic_pairs.append((topic, topic))

# generate list of expected pairwise sub-matrix output files;
# since loading the full matrix is slow, these will be generated in a single rule
pairwise_submats = [join(out_dir, "submats", f"{x[0]}_{x[1]}_comat.mtx") for x in topic_pairs]
pairwise_vocabs = [join(out_dir, "submats", f"{x[0]}_{x[1]}_vocab.feather") for x in topic_pairs]

#
# Snakemake rules
#
# rule all:
#     input:
#         expand(join(out_dir, "submats", "{topic1}_{topic2}.mtx"), topic1=topics, topic2=topics)
        # join(out_dir, "submats", "{topic1}_{topic2}.mtx")

rule build_topic_pair_submatrices:
    input:
        comat=join(out_dir, 'bioconcepts2pubtatorcentral_comat_combined.mtx'),
        vocab=join(out_dir, 'bioconcepts2pubtatorcentral_vocab.feather')
    output:
        submats=pairwise_submats,
        vocabs=pairwise_vocabs
    run:
        # load data
        mat = mmread(input['comat'])
        vocab = pd.read_feather(input['vocab'])
        
        # convert to csr format for faster slicing (better than coo/lil for large
        # dimensions)
        mat = mat.tocsr()

        # iterate over topic pairs and extract relevant sub-matrices for each
        for i in range(len(output['submats'])):
            print("Processing sub-matrix pair %d/%d...", i, len(output['submats']))

            submat_outfile = output['submats'][i]
            vocab_outfile = output['vocabs'][i]

            # extract topics ("topic1_topic2_comat.mtx")
            topic1, topic2, suffix = os.path.basename(submat_outfile).split("_")

            # get indices for each topic
            mask1 = vocab.type == topic1
            mask2 = vocab.type == topic2

            ind1 = vocab.ind[mask1].values
            ind2 = vocab.ind[mask2].values

            # extract relevant rows and columns; must be done one at a time for
            # the sparses matrices
            # submat = mat[ind1].T[ind2].T
            submat = mat[ind1, :][:, ind2]

            mat = mat.tocsr()
            submat = mat[ind1, :][:, ind2]

            # create modified vocab dataframes with updated indices
            vocab_mod = pd.DataFrame({
                "type": list(vocab.type[mask1].values) + list(vocab.type[mask2].values),
                "concept_id": list(vocab.concept_id[mask1].values) + list(vocab.concept_id[mask2].values),
                "ind": list(ind1) + list(ind2)
            })

            # store results
            mmwrite(submat_outfile.replace(".mtx", ""), submat)
            vocab_mod.to_feather(vocab_outfile)

# rule gene_submats:
#     input:
#         comat=join(out_dir, 'bioconcepts2pubtatorcentral_comat_combined.mtx'),
#         vocab=join(out_dir, 'bioconcepts2pubtatorcentral_vocab.feather')
#     output:
#         join(out_dir, 'bioconcepts2pubtatorcentral_gene_submat.mtx')
#     run:
#         # load gene annotations
#         genes = pd.read_csv(config["annotables"], sep="\t")
#
#         # drop entries without entrez ids
#         genes = genes[np.logical_not(np.isnan(genes.entrez))]
#
#         # fix type for entrez column
#         genes['entrez'] = genes.entrez.astype(int).astype(str)
#
#         # load co-occurence matrix and vocab
#         mat = mmread(input['comat'])
#         vocab = pd.read_feather(input['vocab'])
#
#         # convert to lil matrix for faster indexing
#         mat = mat.tolil()
#
#         # store submat-specific indices
#
#
#         # mtor (human mtor entrez id = 2475)
#         # ind = 94299
#
#         # look-up indices associated with concepts of interest
#         concept_ind = []
#
#         breakpoint()
#
#         # to find concept id associated with concept of interest, one do something
#         # like the following with the filtered pubtator dataset:
#         # dat[dat.mentions.str.lower().str.contains('myeloma', regex=False).fillna(False)]
#
#
#         # for concept, concept_id in config['concepts_of_interest'].items():
#         #     myeloma = dat[dat.mentions.str.lower().str.contains('myeloma', regex=False).fillna(False)]
#         #  myeloma.concept_id.value_counts()
#         # ind = 289534
#         #vocab[vocab.concept_id == 'MESH:D009101'].ind.values[0]
#
#         myeloma_ind = 289534
#
#         # iterate over genes and retrieve gene-myeloma counts
#         res = []
#
#         for i in range(genes.shape[0]):
#             entrez_id = genes.iloc[i].entrez
#             symbol = genes.iloc[i].symbol
#
#             gene_matches = vocab.ind[vocab.concept_id == entrez_id]
#
#             # if gene not found in co-matrix, skip
#             if len(gene_matches) == 0:
#                 continue
#
#             gene_ind = gene_matches.values[0]
#
#             res.append({
#                 'gene': symbol,
#                 'count': mat[gene_ind, myeloma_ind]
#             })
#
#         df = pd.DataFrame(res)
#
#         # in cases where multiple entrez ids map to the same gene, choose the max count observed
#         df = df.groupby('gene').max('count')
#
#         df.sort_values('count', ascending=False)
#
#         # BETTER:
#         # get a list of all entrez gene id indices
#         # get a list of all target interaction ids (myeloma, drug sensitivity, biomarker,
#         # etc.)
#         # then select those sets of rows/cols
#
#         # save result
#
#

rule combine_matrices:
    input:
        comatrices=expand(join(out_dir, 'co-occurrence', 'bioconcepts2pubtatorcentral_comat_{i}.npz'),
                          i=range(config['clustering']['num_clusters'])),
        vocab=join(out_dir, 'bioconcepts2pubtatorcentral_vocab.feather')
    output:
        join(out_dir, 'bioconcepts2pubtatorcentral_comat_combined.mtx'),
        join(out_dir, 'bioconcepts2pubtatorcentral_vocab_counts.feather'),
    run:
        mat = sparse.load_npz(input['comatrices'][0])

        for infile in input['comatrices'][1:]:
            mat = mat + sparse.load_npz(infile)

        # determine total counts for each concept id and add to vocab dataframe
        vocab['count'] = mat.sum(axis=0).tolist()[0]
        vocab.to_feather(output[1])

        # make symmetric
        # https://stackoverflow.com/questions/19311353/element-wise-maximum-of-two-sparse-matrices
        mask = (mat > mat.T).astype(int)
        mat = mask.multiply(mat - mat.T) + mat.T

        # store in Matrix Market format for R compatibility
        mmwrite(output[0].replace(".mtx", ""), mat, symmetry='symmetric')

rule filter_dataset:
    input:
        config['pubtator_data_path']
    output:
        join(out_dir, 'bioconcepts2pubtatorcentral_filtered.feather'),
        join(out_dir, 'bioconcepts2pubtatorcentral_vocab.feather')
    run:
        # load data;
        # note: a very small number of malformed-lines containing >4 tab separators
        # were found when processing the data and are explicitly ignored below
        dat = pd.read_csv(input[0], sep='\t', error_bad_lines=False,
                          names=['pmid', 'type', 'concept_id', 'mentions', 'resource'],
                          dtype={'pmid': np.int32, 'type': str, 'concept_id': str,
                                 'mentions': str, 'resource': str})

        pmid_counts = dat.pmid.value_counts()

        # filter out any pubmed articles with a greater than expected number of concept
        # ids assigned to them
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

        # limit to concepts of interest
        dat = dat[dat.type.isin(config['include_concepts'])]

        # TODO: print removal stats
        print(f"Final dataset size: {dat.shape[0]} rows")

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
        vocab_df = dat[['type', 'concept_id', 'ind']].drop_duplicates().sort_values('ind')

        vocab_df.reset_index().drop('index', axis=1).to_feather(output[1])

rule cluster_articles:
    input:
        join(out_dir, 'bioconcepts2pubtatorcentral_filtered.feather')
    output:
        join(out_dir, 'bioconcepts2pubtatorcentral_article_matrix.feather'),
        join(out_dir, 'bioconcepts2pubtatorcentral_article_clusters.feather')
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

        # scale (optional)
        if config['clustering']['scale']:
            X = scale(X)

        # pca projection (optional)
        if config['clustering']['pca']:
            from sklearn.decomposition import PCA

            pca = PCA(n_components=10, random_state=0)
            X = pca.fit_transform(X)

        # NOTE: for the ~13,000,000 x 160 matrix, the sparsity is ~98%.
        kmeans = MiniBatchKMeans(n_clusters=config['clustering']['num_clusters'], random_state=0, batch_size=1000000, max_iter=10).fit(X)

        # store cluster labels
        clusters = pd.DataFrame({'cluster': kmeans.labels_}, index=pmids)

        clusters.reset_index().to_feather(output[1])

rule split_dataset:
    input:
        join(out_dir, 'bioconcepts2pubtatorcentral_filtered.feather'),
        join(out_dir, 'bioconcepts2pubtatorcentral_article_clusters.feather')
    output: expand(join(out_dir, 'subsets', 'bioconcepts2pubtatorcentral_filtered_{i}.feather'), i=range(config['clustering']['num_clusters']))
    run:
        # load full dataset
        dat = pd.read_feather(input[0]).drop(['index'], axis=1)

        pmids = list(dat.pmid.unique())

        # load paper cluster assignments
        clusters = pd.read_feather(input[1]).set_index('pmid')

        # split into n chunks with similar numbers of articles in each
        for i in clusters.cluster.unique():
            # get articles in cluster and store result
            subset_ids = clusters[clusters.cluster == i].index

            # cluster indices start from 0 / snakemake starts from 1
            dat[dat.pmid.isin(subset_ids)].reset_index().to_feather(output[i])

rule build_cooccurrence_matrices:
    input:
        join(out_dir, 'subsets', 'bioconcepts2pubtatorcentral_filtered_{i}.feather'),
        join(out_dir, 'bioconcepts2pubtatorcentral_vocab.feather')
    output:
        join(out_dir, 'co-occurrence', 'bioconcepts2pubtatorcentral_comat_{i}.npz')
    run:
        # load data subset
        dat = pd.read_feather(input[0]).drop(['index'], axis=1)

        # load complete vocabulary
        vocab_df = pd.read_feather(input[1])

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
            if counter % 1000 == 0:
                print("Processing article %d / %d (%0.2f%%)..." % (counter + 1,
                        num_articles, 100 * (counter + 1) / num_articles))

            annots = dat[dat.pmid == pmid].ind.unique()

            # for i, j in itertools.combinations(range(len(annots)), 2):
            for i, j in itertools.combinations(annots, 2):
                mat[i, j] += 1

        t1 = time.time()
        tdelt = (t1 - t0) / (60 * 60)
        rate = num_articles / (t1 - t0)

        print("Finished building co-occurrence matrix.. (Time: %0.1f hours, Rate: %0.2f articles/second)" % (tdelt, rate))

        # convert to coo format and store result
        sparse.save_npz(output[0], mat.tocoo())

