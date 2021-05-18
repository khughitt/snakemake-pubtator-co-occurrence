PubTator Co-occurence Matrix Snakemake Pipeline
===============================================

V. Keith Hughitt (Sept 2020)

Overview
--------

This pipeline serves to generate a concept x concept co-occurence matrix using the bulk
data downloads available from the [PubTator Central FTP site](ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral).

The [Snakemake](https://snakemake.readthedocs.io/en/stable/) computational pipeline is
used, enabling the pipeline to be easily deployed both locally, as well as on a
high-performancing computing environment such as a cluster.

The basic steps in the pipeline are:

1. Filter bioconcepts to remove unneeded and low-frequency concepts
2. Cluster pubmed articles based on their usage of a subset of high-frequency concepts
3. Divide dataset into N subsets based on the clusters determined above
4. Compute co-occurence matrices separately for each data subset
5. Combine sub co-occurrence matrices into a final complete matrix

Usage
-----

Create and activate a conda environment with the neccessary components:

```
conda create -n pubtator --file requirements.txt
conda activate pubtator
```

Download the most recent version of the full abstract-level PubTator bioconcept
annotations file (`bioconcepts2pubtatorcentral.gz`) from: 

[ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral](ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral)

Copy and modify the example configuration file located in the `config/` directory,
specifying the location of the PubTator data data downloaded above and modifying the
pipeline parameters as desired.

To run the pipeline with a total of four threads, one can then do:

```
snakemake -j 4 --config-file path/to/config.yml
```

Pubtator Data Statistics
------------------------

Some useful statistics about that annotations present in full PubTator bioconcepts dump.

These numbers are based on the most recent version of `bioconcepts2pubtatorcentral.gz`,
at the time of writing, August 13, 2020.

Unique values present in each column:

```
pmid          23148837
type                 9
concept_id     9180086
mentions      25818648
resource           955
```

Frequencies of each concept "type":

```
Disease        96979657
Chemical       89112532
Gene           50087304
Species        36151864
Mutation        2746529
CellLine        1939887
Genus              3934
Strain              919
DomainMotif         429
```

Distribution of low-frequency concept ids:

```
freq  count
1     7525737
2      757479
3      269939
4      109938
5       73627
6       46751
7       38308
8       28303
9       21778
10      21298
```

Based on the above, we can see that there are 7,525,737 annotations which only appear
once across all annotated articles.

Notes
-----

### Identifiers

The same mention can be mapped to multiple concept ids:

```
               pmid  type concept_id mentions   resource   ind
1118       27028000  Gene      56717     mTOR  GNormPlus   827
3732       31247000  Gene       2475     mTOR  GNormPlus  2291
...
```

In the above, '56717' refers to mouse mTOR while 2475 refers to human mTOR:

- https://www.ncbi.nlm.nih.gov/gene/56717
- https://www.ncbi.nlm.nih.gov/gene/2475

Similarly same concept id can be used to refer to different mentions, e.g.:

```
dat[dat.concept_id == '2475'].groupby('mentions').concept_id.value_counts().sort_values(ascending=False)

mentions                                                               concept_id
mTOR                                                                   2475          38573
mammalian target of rapamycin|mTOR                                     2475          19431
mammalian target of rapamycin                                          2475           3821
FRAP                                                                   2475           1777
...
```
