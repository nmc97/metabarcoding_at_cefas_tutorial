=============
Metabarcoding
=============

What is metabarcoding?
======================

# insert description and/or links

Process:
========

1. Set up experiment - keep the rest of the analysis in mind when setting up the experiment

  - choose a target sequence

2. sequencing

  - note: Nanopore metabarcoding while produces longer sequences, the databases for assigning taxonomy to these sequences may be poor in comparison to other metabarcoding approaches

3. Data archiving


Quality Control: Read quality
=============================

Forward and reverse reads in MiSeq - can struggle on metabarcoding runs - low sequence diversity
does better on forward than reverse reads - lower quality
Decide if you should use reverse reads at all

Overlapping reads - do they overlap enough? If not how do I include them


OTU's ASV's
===========

Before deciding what clustering method to use it is important tto understand the different types of clusters that you may want to produce.

[insert]

Clustering
==========

There are three types of algorithm for clustering metabarcoding reads into OTU's or ASV's:

`De novo clustering - threshold <2\. De novo Clustering strategy - defined threshold_>`_

1\. Alignment-based strategy
----------------------------

Aligns reads to a database

**Software:**

* `Kraken <https://github.com/DerrickWood/kraken2/wiki/Manual>`_
* `Centrifuge <http://ccb.jhu.edu/software/centrifuge/>`_
* `Minimap <https://github.com/lh3/minimap2>`_

Option: Visualise with `Pavian <https://github.com/fbreitwieser/pavian>`_

**Considerations:**

* Can miss uncharacterised species
* minimap alignments require filtering (parameters are questionable)

2\. De novo Clustering strategy - defined threshold
---------------------------------------------------

Picks a threshold at which to define a cluster - not really taxa arbitrary grouping

Considerations:
Been in use for a long time so well understood - used for a long time
More parameters than alignment strategy so more complicated

3\. Clustering with guided clustering instead of thresholds
-----------------------------------------------------------

Not as arbitrary as threshold-based analysis
Accounts for sequencing errors

* `USEARCH <http://www.drive5.com/usearch/>`_
* `SWARM <https://github.com/torognes/swarm>`_
* `DADA2 <https://benjjneb.github.io/dada2/>`_

Identify Chimeras
=================

1.  What proportion foo the reads align to the reference?
2.  Chimera could be 2 species you haven't seen before
3.  Check OTU's individually
4.  Check against reference
5.  More abundant OTU's more likely to be real
6.  OTU's for every library and them split

R

Useful databases
=================

Choosing a database:

- SILVA
- PR2 - `18S database <https://pr2-database.org/>`_

Cox-1 gene databases:

- `Automated high throughput animal CO1 metabarcoding classification <https://www.nature.com/articles/s41598-018-22505-4>`_

Statistics
==========

Don't do this on POD
Phyloseq is good but is limited due to the developer
Output files/ abundance file - try to have them in `.biom` format - relatively universal

Alpha diversity
^^^^^^^^^^^

Measure of abundance of species

Lots of metrics!

- Count number of Taxa
- Treat as a sample of the overal population and attempt to calculate the population - chao
- level of evenness - how evenly they split

Normalising
^^^^^^^^^^^

- accounting for sequencing depth before doing your analysis
- subsample seq dataset
    - check multiple coverage levels and plot to see if it levels off
    - more sequences - more errors so more OTU's line will never be flat
    - accounts for sequencing depth twice

Beta diversity
^^^^^^^^^^^^^^

- unsupervised analysis (doesn't know which samples are in which group)
- based on the abundance of these ...
- Do these cluster together or apart
- maximised variability
- normal to do a lot of normalising before this step
    - lots of different ways


Packages for statistics
^^^^^^^^^^^^^^^^^^^^^^^

Process - Dada2
trimming, deduplicating, clustering
chimeras
Taxon Id
Table of Abundances
