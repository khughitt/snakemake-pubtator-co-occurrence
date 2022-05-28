"""
Pubtator Co-occurrence matrix pipeline (v2)

KH (May 18, 2021)

Processes PubtatorCentral annotations and generates some useful "clean" versions of the
data, along with some co-occurrence matrices relating to different concepts, for Human.
"""
import itertools
import json
import os
import random
import time
import numpy as np
import pandas as pd
from os.path import join
from scipy import sparse
from scipy.io import mmread, mmwrite
from sklearn.preprocessing import scale

configfile: "config/config.yml"

random.seed(config['random_seed'])

out_dir = join(config["out_dir"], config["version"])

if config['debug']:
    out_dir = join(out_dir, 'sampled')

# encoder to convert int64 elements to generic ints and sets to lists during
# json serialization
def encoder(object):
    if isinstance(object, np.generic):
        return object.item()
    elif isinstance(object, set):
        return list(object)

#
# Snakemake rules
#
rule all:
    input:
        expand(join(out_dir, 'all-genes/chemicals/{chemical}.feather'), chemical=config['chemicals'].keys()),
        expand(join(out_dir, 'all-chemicals/genes/{gene}.feather'), gene=config['genes'].keys()),
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_diseases.feather'),
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_chemicals.feather'),
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_genes.feather'),
        join(out_dir, 'pmids/gene-pmids.json'),
        join(out_dir, 'pmids/chemical-pmids.json'),
        join(out_dir, 'pmids/disease-pmids.json'),
        join(out_dir, 'co-occurrence/genes.feather'),
        join(out_dir, 'co-occurrence/diseases.feather'),
        join(out_dir, 'co-occurrence/chemicals.feather'),
        join(out_dir, 'co-occurrence/genes-diseases.feather'),
        join(out_dir, 'co-occurrence/genes-chemicals.feather'),
        join(out_dir, 'identifiers/human_entrez_ids.txt'),
        join(out_dir, 'gene-counts/pmid_gene_counts.feather'),
        join(out_dir, 'metadata/mesh_disease_mapping.tsv'),
        join(out_dir, 'metadata/mesh_chemical_mapping.tsv'),
        # join(out_dir, 'pmids/pmid-genes.json'),

rule gene_disease_comat:
    input:
        join(out_dir, 'pmids/gene-pmids.json'),
        join(out_dir, 'pmids/disease-pmids.json')
    output:
        join(out_dir, 'co-occurrence/genes-diseases.feather')
    run:
        # load gene and drug pmid mappings
        with open(input[0]) as fp:
            all_gene_pmids = json.load(fp)

        with open(input[1]) as fp:
            all_disease_pmids = json.load(fp)

        # create empty matrix to store gene-disease co-occurrence counts
        entrez_ids = all_gene_pmids.keys()
        num_genes = len(entrez_ids)

        mesh_ids = all_disease_pmids.keys()
        num_diseases = len(mesh_ids)

        comat = np.empty((num_genes, num_diseases))
        comat.fill(np.nan)

        # iterate over pairs of genes
        for i, gene in enumerate(entrez_ids):
            # get pubmed ids associated with gene
            gene_pmids = all_gene_pmids[gene]

            for j, disease in enumerate(mesh_ids):
                # get pubmed ids associated with disease
                disease_pmids = all_disease_pmids[disease]

                # compute gene-disease co-occurrence count
                num_shared = len(set(gene_pmids).intersection(disease_pmids))

                comat[i, j] = num_shared

        # store gene-disease co-occurrence matrix
        comat = pd.DataFrame(comat, index=entrez_ids, columns=mesh_ids)
        comat.reset_index().rename(columns={'index': 'entrez_id'}).to_feather(output[0])

rule gene_chemical_comat:
    input:
        join(out_dir, 'pmids/gene-pmids.json'),
        join(out_dir, 'pmids/chemical-pmids.json')
    output:
        join(out_dir, 'co-occurrence/genes-chemicals.feather')
    run:
        # load gene and drug pmid mappings
        with open(input[0]) as fp:
            all_gene_pmids = json.load(fp)

        with open(input[1]) as fp:
            all_chemical_pmids = json.load(fp)

        # create empty matrix to store gene-chemical co-occurrence counts
        entrez_ids = all_gene_pmids.keys()
        num_genes = len(entrez_ids)

        mesh_ids = all_chemical_pmids.keys()

        # exclude "-" and chemicals with fewer associated pmids
        # TODO: move upstream prior to next major run..

        # note: in current version of pipeline, setting a minimum pmid count of "50" for
        # the chemical data results in 17,811/77,041 chemicals being included.. 
        chem_pmid_counts = pd.Series([len(a) for a in all_chemical_pmids.values()],
                                     index=all_chemical_pmids.keys())

        MIN_PMID_COUNT = 50

        to_keep = chem_pmid_counts[chem_pmid_counts >= MIN_PMID_COUNT].index
        mesh_ids = [x for x in mesh_ids if x in to_keep]

        mesh_ids.remove('-')

        num_chemicals = len(mesh_ids)

        comat = np.empty((num_genes, num_chemicals))
        comat.fill(np.nan)

        print(f"Processing {num_genes} genes...")

        # iterate over pairs of genes
        for i, gene in enumerate(entrez_ids):
            # get pubmed ids associated with gene
            gene_pmids = all_gene_pmids[gene]

            if i % 100 == 0:
                print(f"gene {i}/{num_genes}...")

            for j, chemical in enumerate(mesh_ids):
                # get pubmed ids associated with chemical
                chemical_pmids = all_chemical_pmids[chemical]

                # compute gene-chemical co-occurrence count
                num_shared = len(set(gene_pmids).intersection(chemical_pmids))

                comat[i, j] = num_shared

        # store gene-chemical co-occurrence matrix
        comat = pd.DataFrame(comat, index=entrez_ids, columns=mesh_ids)
        comat.reset_index().rename(columns={'index': 'entrez_id'}).to_feather(output[0])

rule gene_gene_comat:
    input:
        join(out_dir, 'pmids/gene-pmids.json')
    output:
        join(out_dir, 'co-occurrence/genes.feather')
    run:
        # load gene/pmid mapping
        with open(input[0]) as fp:
            gene_pmids = json.load(fp)

        # create empty matrix to store gene-gene co-occurrence counts
        entrez_ids = gene_pmids.keys()
        num_genes = len(entrez_ids)

        comat = np.empty((num_genes, num_genes))
        comat.fill(np.nan)

        # iterate over pairs of genes
        for i, gene1 in enumerate(entrez_ids):
            # get pubmed ids associated with gene 1
            gene1_pmids = gene_pmids[gene1]

            for j, gene2 in enumerate(entrez_ids):
                # skip symmetric comparisons
                if not np.isnan(comat[i, j]):
                    continue

                gene2_pmids = gene_pmids[gene2]

                # compute gene-gene co-occurrence count
                num_shared = len(set(gene1_pmids).intersection(gene2_pmids))

                comat[i, j] = comat[j, i] = num_shared

        # store gene-gene co-occurrence matrix
        comat = pd.DataFrame(comat, index=entrez_ids, columns=entrez_ids)
        comat.reset_index().rename(columns={'index': 'entrez_id'}).to_feather(output[0])

rule disease_disease_comat:
    input:
        join(out_dir, 'pmids/disease-pmids.json')
    output:
        join(out_dir, 'co-occurrence/diseases.feather')
    run:
        # load disease/pmid mapping
        with open(input[0]) as fp:
            disease_pmids = json.load(fp)

        # create empty matrix to store disease-disease co-occurrence counts
        mesh_ids = disease_pmids.keys()
        num_diseases = len(mesh_ids)

        comat = np.empty((num_diseases, num_diseases))
        comat.fill(np.nan)

        # iterate over pairs of diseases
        for i, disease1 in enumerate(mesh_ids):
            # get pubmed ids associated with disease 1
            disease1_pmids = disease_pmids[disease1]

            for j, disease2 in enumerate(mesh_ids):
                # skip symmetric comparisons
                if not np.isnan(comat[i, j]):
                    continue

                disease2_pmids = disease_pmids[disease2]

                # compute disease-disease co-occurrence count
                num_shared = len(set(disease1_pmids).intersection(disease2_pmids))

                comat[i, j] = comat[j, i] = num_shared

        # store disease-disease co-occurrence matrix
        comat = pd.DataFrame(comat, index=mesh_ids, columns=mesh_ids)
        comat.reset_index().rename(columns={'index': 'mesh_id'}).to_feather(output[0])

rule disease_pmid_mapping:
    input:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_diseases.feather')
    output:
        join(out_dir, 'pmids/disease-pmids.json')
    run:
        # load disease data
        disease_dat = pd.read_feather(input[0])

        # iterate over diseases
        mesh_ids = list(disease_dat.concept_id.unique())
        num_diseases = len(mesh_ids)

        # iterate over diseases and retrieve associated pubmed ids for each
        disease_pmids = {}

        for mesh_id in mesh_ids:
            mask = disease_dat.concept_id == mesh_id
            disease_pmids[mesh_id] = set(disease_dat[mask].pmid.values)

        # store disease -> pmid mapping as json
        with open(output[0], "w") as fp:
            fp.write(json.dumps(disease_pmids, default=encoder))

rule chemical_chemical_comat:
    input:
        join(out_dir, 'pmids/chemical-pmids.json')
    output:
        join(out_dir, 'co-occurrence/chemicals.feather')
    run:
        # load chemical/pmid mapping
        with open(input[0]) as fp:
            unfiltered_chemical_pmids = json.load(fp)

        # filter any chemicals with less than N citations; helps keep the co-occurrence
        # matrix size more manageable while only losing relatively low-information
        # low-citation counts..
        chemical_pmids = {}

        for mesh_id in unfiltered_chemical_pmids:
            if len(unfiltered_chemical_pmids[mesh_id]) > config['concept_id_min_freq_chem']:
                chemical_pmids[mesh_id] = unfiltered_chemical_pmids[mesh_id]

        # create empty matrix to store chemical-chemical co-occurrence counts
        mesh_ids = chemical_pmids.keys()
        num_chemicals = len(mesh_ids)

        comat = np.empty((num_chemicals, num_chemicals))
        comat.fill(np.nan)

        # iterate over pairs of chemicals
        for i, chemical1 in enumerate(mesh_ids):
            # get pubmed ids associated with chemical 1
            chemical1_pmids = chemical_pmids[chemical1]

            for j, chemical2 in enumerate(mesh_ids):
                # skip symmetric comparisons
                if not np.isnan(comat[i, j]):
                    continue

                chemical2_pmids = chemical_pmids[chemical2]

                # compute chemical-chemical co-occurrence count
                num_shared = len(set(chemical1_pmids).intersection(chemical2_pmids))

                comat[i, j] = comat[j, i] = num_shared

        # store chemical-chemical co-occurrence matrix
        comat = pd.DataFrame(comat, index=mesh_ids, columns=mesh_ids)
        comat.reset_index().rename(columns={'index': 'mesh_id'}).to_feather(output[0])

rule chemical_pmid_mapping:
    input:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_chemicals.feather')
    output:
        join(out_dir, 'pmids/chemical-pmids.json')
    run:
        # load chemical data
        chemical_dat = pd.read_feather(input[0])

        # iterate over chemicals
        mesh_ids = list(chemical_dat.concept_id.unique())
        num_chemicals = len(mesh_ids)

        # iterate over chemicals and retrieve associated pubmed ids for each
        chemical_pmids = {}

        for mesh_id in mesh_ids:
            mask = chemical_dat.concept_id == mesh_id
            chemical_pmids[mesh_id] = set(chemical_dat[mask].pmid.values)

        # store chemical -> pmid mapping as json
        with open(output[0], "w") as fp:
            fp.write(json.dumps(chemical_pmids, default=encoder))

rule gene_pmid_mapping:
    input:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_genes.feather')
    output:
        join(out_dir, 'pmids/gene-pmids.json')
    run:
        # load gene data
        gene_dat = pd.read_feather(input[0])

        # iterate over genes
        entrez_ids = list(gene_dat.concept_id.unique())
        num_genes = len(entrez_ids)

        # iterate over genes and retrieve associated pubmed ids for each
        gene_pmids = {}

        for entrez_id in entrez_ids:
            mask = gene_dat.concept_id == entrez_id
            gene_pmids[entrez_id] = set(gene_dat[mask].pmid.values)

        # store gene -> pmid mapping as json
        with open(output[0], "w") as fp:
            fp.write(json.dumps(gene_pmids, default=encoder))

#
# create a mapping from pmid to genes;
# this can be used to exclude or down-weight articles which reference large numbers of
# genes, and thus, and not very specific / informative.
#
# rule pmid_gene_mapping:
#     input:
#         join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_genes.feather')
#     output:
#         join(out_dir, 'pmids/pmid-genes.json')
#     run:
#         # load gene data
#         gene_dat = pd.read_feather(input[0])
#
#         # iterate over pubmed ids
#         pmids = list(gene_dat.pmid.unique())
#         num_pmids = len(pmids)
#
#         # iterate over pubmed articles, and get associated genes.
#         pmid_genes = {}
#
#         # for entrez_id in entrez_ids:
#         for pmid in pmids:
#             mask = gene_dat.pmid == pmid
#             pmid_genes[pmid] = set(gene_dat[mask].concept_id.values)
#
#         # store gene -> pmid mapping as json
#         with open(output[0], "w") as fp:
#             fp.write(json.dumps(pmid_genes, default=encoder))

#
# computes chemical-gene co-occurrence counts for:
#
# - all chemicals
# - specific genes
#
rule chemicals_to_genes:
    input:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_genes.feather'),
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_chemicals.feather')
    output:
        join(out_dir, 'all-chemicals/genes/{gene}.feather')
    params:
        gene='{gene}'
    run:
        # load gene data
        gene_dat = pd.read_feather(input[0])

        # get pubmed for articles referencing the gene
        gene_entrez = config['genes'][params.gene]
        mask = gene_dat.concept_id == gene_entrez

        gene_pmids = set(gene_dat[mask].pmid.values)

        # load chemical data
        chemical_dat = pd.read_feather(input[1])

        # exclude chemicals with no associated mesh id
        chemical_dat = chemical_dat.loc[chemical_dat.concept_id != '-']

        # limit to pubtator entries for articles referencing gene
        chemical_dat = chemical_dat[chemical_dat.pmid.isin(gene_pmids)]

        # get counts for each chemical
        counts = chemical_dat[['concept_id']].value_counts().sort_values(ascending=False)

        # clean-up dataframe and add chemical symbol, etc.
        counts = counts.reset_index()
        counts.columns = ['mesh_id', 'n']

        # load mesh/chemical mapping
        mapping = pd.read_csv("data/mesh_chemical_mapping.tsv", sep="\t")
        mapping.columns = ['mesh_id', 'chemical']

        # for multi-mapped entrez chemicals, simply choose the first match from each group 
        # as a representative
        counts = counts.merge(mapping, on=['mesh_id'])

        # store result
        counts.to_feather(output[0])

#
# computes chemical-gene co-occurrence counts for:
#
# - all genes
# - specific chemicals
#
rule genes_to_chemicals:
    input:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_chemicals.feather'),
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_genes.feather')
    output:
        join(out_dir, 'all-genes/chemicals/{chemical}.feather')
    params:
        chemical='{chemical}'
    run:
        # load chemical data
        chemical_dat = pd.read_feather(input[0])

        # get pubmed for articles referencing the chemical
        mesh_id = config['chemicals'][params.chemical]
        mask = chemical_dat.concept_id == mesh_id

        pmids = set(chemical_dat[mask].pmid.values)

        # load gene data
        gene_dat = pd.read_feather(input[1])

        # limit to entries for articles referencing chemical
        gene_dat = gene_dat[gene_dat.pmid.isin(pmids)]

        # convert from long to wide
        #gene_dat['present'] = 1
        #dat_wide = gene_dat.pivot_table(index='concept_id', columns='pmid', values='present', fill_value=0)

        # get counts for each gene
        counts = gene_dat[['concept_id']].value_counts().sort_values(ascending=False)

        # clean-up dataframe and add gene symbol, etc.
        counts = counts.reset_index()
        counts.columns = ['entrez', 'n']

        # load gene annotations
        grch37 = pd.read_csv("data/entrez_gene_mapping.tsv", sep="\t")

        # drop genes with missing entrez ids and fix type for remaining ids
        grch37 = grch37[~grch37.entrez.isna()]
        grch37 = grch37.astype({'entrez': int}).astype({'entrez': str})

        # for multi-mapped entrez genes, simply choose the first match from each group 
        # as a representative
        grch37 = grch37.groupby('entrez').first()

        counts = counts.merge(grch37, on=['entrez'])

        # store result
        counts.to_feather(output[0])

# count number of unique genes associated with each pmid, after filtering..
rule pubmed_gene_counts:
    input:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_genes.feather')
    output:
        join(out_dir, 'gene-counts/pmid_gene_counts.feather')
    run:
        dat = pd.read_feather(input[0])
 
        counts = dat.groupby('pmid').concept_id.nunique()
        counts.rename("num").reset_index().to_feather(output[0])

rule create_disease_subset:
    input:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human.tsv')
    output:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_diseases.feather')
    run:
        # load filtered dataset
        dat = pd.read_csv(input[0], sep='\t')

        # create and store subset with only disease-related entries
        dat = dat[(dat.type == 'Disease')]
        dat.reset_index(drop=True).to_feather(output[0])

rule create_chemical_subset:
    input:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human.tsv')
    output:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_chemicals.feather')
    run:
        # load filtered dataset
        dat = pd.read_csv(input[0], sep='\t')

        # create and store subset with only chemical-related entries
        dat = dat[(dat.type == 'Chemical')]
        dat.reset_index(drop=True).to_feather(output[0])

rule create_gene_subset:
    input:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human.tsv'),
        join(out_dir, 'identifiers/human_entrez_ids.txt')
    output:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_genes.feather')
    run:
        # load filtered dataset
        dat = pd.read_csv(input[0], sep='\t')

        # load list of human entrez gene identifiers
        with open(input[1]) as fp:
            entrez_ids = [x.strip() for x in fp.readlines()]

        # create a subset with only gene entries, restricted to GNormPlus/Human Entrez
        # gene identifiers..
        entrez_mask = dat.concept_id.isin(entrez_ids)

        dat = dat[(dat.type == 'Gene') & (dat.resource == 'GNormPlus') & entrez_mask]

        # store gene subset
        dat.reset_index(drop=True).to_feather(output[0])

#
# note: currently saving filtered data as uncompressed .tsv file instead of feather due
# to the higher memory usage for writing compressed feather files. (kh may21)
#
rule filter_dataset:
    input:
        config['pubtator_data_path']
    output:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human.tsv')
    run:
        # load data
        dat = pd.read_csv(input[0], sep='\t',
                          names=['pmid', 'type', 'concept_id', 'mentions', 'resource'],
                          dtype={'pmid': np.int32, 'type': str, 'concept_id': str,
                                 'mentions': str, 'resource': str})

        # if debug mode is enabled, sub-sample the data to speed up development
        if config['debug']:
            print(f"[DEBUG] Randomly sampling {config['sample_pmids']} pubmed articles...")

            sampled_pmids = np.random.choice(dat.pmid.unique(), config['sample_pmids'])
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

        num_kept = 0 if True not in mask_counts else mask_counts[True]
        num_total = dat.shape[0]

        num_dropped = 0 if False not in mask_counts else mask_counts[False]

        pct_dropped = 100 * (num_dropped / num_total)

        print("Excluding %d / %d (%0.2f%%) concept ids which appear <%d times.." % (num_dropped, num_total, pct_dropped, config['concept_id_min_freq']))

        dat = dat[mask]

        # print summary of filtering results
        print(f"Final dataset size: {dat.shape[0]} rows ({len(set(dat.pmid))} articles)")
        print(f"Removed {num_entries_orig} rows in total ({num_pmids_orig} articles)") 

        # save filtered dataset
        #dat.reset_index(drop=True).to_csv(output[0], sep='\t')
        dat.to_csv(output[0], index=False, sep='\t')

rule create_mesh_chem_mapping:
    input:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_chemicals.feather')
    output:
        join(out_dir, 'metadata/mesh_chemical_mapping.tsv'),
    params:
        field='chemical'
    script:
        "scripts/create_mesh_mapping.R"

rule create_mesh_disease_mapping:
    input:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_diseases.feather')
    output:
        join(out_dir, 'metadata/mesh_disease_mapping.tsv'),
    params:
        field='disease'
    script:
        "scripts/create_mesh_mapping.R"

rule get_entrez_ids:
    output:
        join(out_dir, 'identifiers/human_entrez_ids.txt')
    script:
        "scripts/create_entrez_human_id_list.R"
