Metabarcoding
=============

What is metabarcoding?
^^^^^^^^^^^^^^^^^^^^^^

# insert description and/or links

General Process:
^^^^^^^^^^^^^^^^

#. Set up experiment - keep the rest of the analysis in mind when setting up the experiment
Choose a target sequence

#. Sequencing

note: Nanopore metabarcoding while produces longer sequences, the databases for assigning taxonomy to these sequences may be poor in comparison to other metabarcoding approaches

#. Data archiving

Read files (fastq.gz, metadata)

#. Quality control

 - trimming - remove low quality reads, adapters and trim low quality read ends.
 - deduplicating - removing reads which appear more than once in the dataset (reads can be overalappinga and highly similar but reads that are identical are redundant)

#. Clustering

  - lots of different ways to do this.
  - What types of clusters are you looking for? OUT's? ASV's?
  - What type of data do you have? e.g: 16S, 18S, COI, long/short reads.

#. Check for chimeras

  - these are sequences which are atrifical

#. Classify Taxa

 - Choose a database based on the organisms and target sequences you are working with

#.  Abundance statistics - alpha and beta Metrics

#. Differential abundance

Quality Control: Read quality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Forward and reverse reads in MiSeq - can struggle on metabarcoding runs - low sequence diversity
does better on forward than reverse reads - lower quality
Decide if you should use reverse reads at all

Overlapping reads - do they overlap enough? If not how do I include them

[insert]

Clustering
^^^^^^^^^^

**OTU's ASV's**

Before deciding what clustering method to use it is important tto understand the different types of clusters that you may want to produce.


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

*3\. Clustering with guided clustering instead of thresholds*
-------------------------------------------------------------
Not as arbitrary as threshold-based analysis
Accounts for sequencing errors

* `USEARCH <http://www.drive5.com/usearch/>`_
* `SWARM <https://github.com/torognes/swarm>`_
* `DADA2 <https://benjjneb.github.io/dada2/>`_

Identify Chimeras
^^^^^^^^^^^^^^^^^

1.  What proportion foo the reads align to the reference?
2.  Chimera could be 2 species you haven't seen before
3.  Check OTU's individually
4.  Check against reference
5.  More abundant OTU's more likely to be real
6.  OTU's for every library and them split

Taxanomic assignment:
^^^^^^^^^^^^^^^^^^^^^

Useful databases
----------------

Choosing a database ...

- SILVA
- PR2 - `18S database <https://pr2-database.org/>`_

Cox-1 gene databases:
- costom database: `DUFA <github.com/uit-metabarcoding/DUFA>`_ : `Paper <https://academic.oup.com/icesjms/article/78/9/3342/6360557#323435484>`_

- `Automated high throughput animal CO1 metabarcoding classification <https://www.nature.com/articles/s41598-018-22505-4>`_

Statistics
^^^^^^^^^^

Don't do this on POD

Phyloseq is good but is limited due to the developer

Output files/ abundance file - try to have them in `.biom` format - relatively universal


Normalising
-----------

- Accounting for sequencing depth before doing your analysis

- Subsample seq dataset

    - Check multiple coverage levels and plot to see if it levels off

    - More sequences - more errors so more OTU's line will never be flat

    - Accounts for sequencing depth twice

Resources:
`Normalization and microbial differential abundance strategies depend upon data characteristics <https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y>`_

-

Alpha diversity
---------------

Alpha diversity is a measure of species abundance in each sample, or all samples pooled.

There a lot's of different metrics which can be used to calculate this; thus, alpha metrics cannot readily be compared between studies.

Metrics:

- Count number of Taxa

- Treat as a sample of the overal population and attempt to calculate the population - Chao

- Level of evenness - how evenly they split

Beta diversity
--------------

- Unsupervised analysis (doesn't know which samples are in which group)

- Based on the abundance of these ...

- Do these cluster together or apart

- Maximised variability

- It is normal to do a lot of normalising before this step

    - lots of different ways to do this

Differential abundance
----------------------

`Microbiome differential abundance methods produce different results across 38 datasets <https://www.nature.com/articles/s41467-022-28034-z>`_

Packages for statistics
-----------------------

Many ststistical packages in r for metabarcoding statistics use the package `vegan <http://vegan.r-forge.r-project.org/FAQ-vegan.html#What-is-vegan_003f>_`
