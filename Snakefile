"""
Pubtator Co-occurrence matrix pipeline (v2)

KH (May 18, 2021)

Processes PubtatorCentral annotations and generates some useful "clean" versions of the
data, along with some co-occurrence matrices relating to different concepts, for Human.
"""
import random
from os.path import join

configfile: "config/config.yml"

random.seed(config['random_seed'])

out_dir = join(config["out_dir"], config["version"])

if config['debug']:
    out_dir = join(out_dir, 'sampled')

#
# Snakemake rules
#
rule all:
    input:
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

rule create_gene_disease_comat:
    input:
        join(out_dir, 'pmids/gene-pmids.json'),
        join(out_dir, 'pmids/disease-pmids.json')
    output:
        join(out_dir, 'co-occurrence/genes-diseases.feather')
    script:
        "scripts/create_gene_disease_comat.py"

rule create_gene_chemical_comat:
    input:
        join(out_dir, 'pmids/gene-pmids.json'),
        join(out_dir, 'pmids/chemical-pmids.json')
    output:
        join(out_dir, 'co-occurrence/genes-chemicals.feather')
    script:
        "scripts/create_gene_chemical_comat.py"

rule create_gene_gene_comat:
    input:
        join(out_dir, 'pmids/gene-pmids.json')
    output:
        join(out_dir, 'co-occurrence/genes.feather')
    script:
        "scripts/create_gene_gene_comat.py"

rule create_disease_disease_comat:
    input:
        join(out_dir, 'pmids/disease-pmids.json')
    output:
        join(out_dir, 'co-occurrence/diseases.feather')
    script:
        "scripts/create_disease_disease_comat.py"

rule create_disease_pmid_mapping:
    input:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_diseases.feather')
    output:
        join(out_dir, 'pmids/disease-pmids.json')
    script:
        "scripts/create_disease_pmid_mapping.py"

rule create_chemical_chemical_comat:
    input:
        join(out_dir, 'pmids/chemical-pmids.json')
    output:
        join(out_dir, 'co-occurrence/chemicals.feather')
    script:
        "scripts/create_chemical_chemical_comat.py"

rule create_chemical_pmid_mapping:
    input:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_chemicals.feather')
    output:
        join(out_dir, 'pmids/chemical-pmids.json')
    script:
        "scripts/create_chemical_pmid_mapping.py"

rule create_gene_pmid_mapping:
    input:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_genes.feather')
    output:
        join(out_dir, 'pmids/gene-pmids.json')
    script:
        "scripts/create_gene_pmid_mapping.py"

rule create_article_gene_counts:
    input:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_genes.feather')
    output:
        join(out_dir, 'gene-counts/pmid_gene_counts.feather')
    script:
        "scripts/create_article_gene_counts.py"

rule create_disease_subset:
    input:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human.tsv')
    output:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_diseases.feather')
    script:
        "scripts/create_disease_subset.py"

rule create_chemical_subset:
    input:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human.tsv')
    output:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_chemicals.feather')
    script:
        "scripts/create_chemical_subset.py"

rule create_gene_subset:
    input:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human.tsv'),
        join(out_dir, 'identifiers/human_entrez_ids.txt')
    output:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human_genes.feather')
    script:
        "scripts/create_gene_subset.py"

# uses externally-generated pubmed embeddings to split articles into batches based on
# article similarity in the embedding space
rule chunk_articles:

# note: currently saving filtered data as uncompressed .tsv file instead of feather due
# to the higher memory usage for writing compressed feather files. (kh may21)
rule filter_bioconcepts:
    input:
        config['pubtator_data_path']
    output:
        join(out_dir, 'filtered/bioconcepts2pubtatorcentral_filtered_human.tsv')
    script:
        "scripts/filter_bioconcepts.py"

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

rule create_entrez_human_id_list:
    output:
        join(out_dir, 'identifiers/human_entrez_ids.txt')
    script:
        "scripts/create_entrez_human_id_list.R"
