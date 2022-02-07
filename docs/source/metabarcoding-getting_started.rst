Metabarcoding
=============

What is metabarcoding?
^^^^^^^^^^^^^^^^^^^^^^

# insert description and/or links

General Process:
^^^^^^^^^^^^^^^^

1. Set up experiment
Keep the rest of the analysis in mind when setting up the experiment
Choose a target sequence

2. Sequencing

note: Nanopore metabarcoding while produces longer sequences, the databases for assigning taxonomy to these sequences may be poor in comparison to other metabarcoding approaches

3. Data archiving

Read files (fastq.gz, metadata)

4. Quality control

- trimming - remove low quality reads, adapters and trim low quality read ends.
- deduplicating - removing reads which appear more than once in the dataset (reads can be overalappinga and highly similar but reads that are identical are redundant)
 e.g.: in Dada2 tutorial

.. code::

  plotQualityProfile(fnFs[1:2]) # plot forward reads
  plotQualityProfile(fnRs[1:2]) # plot reverse reads

  # Place filtered files in filtered/ subdirectory
  filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names

  # decide to truncate forward reads at 240 and reverse reads at 160: truncLen=c(240,160)
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

5. Clustering

- lots of different ways to do this.
- What types of clusters are you looking for? OUT's? ASV's?
- What type of data do you have? e.g.: 16S, 18S, COI, long/short reads.

6. Check for chimeras

- these are sequences which are artificial

7. Classify Taxa

- Choose a database based on the organisms and target sequences you are working with

8. Abundance statistics - alpha and beta Metrics

9. Differential abundance Analysis

Quality Control: Read quality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Forward and reverse reads in MiSeq - can struggle on metabarcoding runs - low sequence diversity
does better on forward than reverse reads - lower quality
Decide if you should use reverse reads at all

Overlapping reads - do they overlap enough? If not how do I include them

[insert]

Programs to use: Fastqc, Trimmomatic, ??? what is standard in metabarcoding?

.. note::

  *Questions to consider:*

  What type of reads do I have?
  Will the forward and reverse reads (if paired) overlap?
  What clustering method will I be using (some account for error in reads somewhat so trimming may be less necessary)
  How many reads do I have?
  What do the read quality checks tell me about the data quality?
  After trimming/ filtering how many reads are there per sample.
  Are there samples that need to be excluded from the rest of the analyses?


Clustering
^^^^^^^^^^

**OTU's vs ASV's**

Before deciding what clustering method to use it is important to understand the different types of clusters that you may want to produce.

*Definitions:*

OTU = Operating Taxonomic Units
ASV = Amplicon Sequence Variant

See: `MICROBIOME INFORMATICS: OTU VS. ASV <https://www.zymoresearch.com/blogs/blog/microbiome-informatics-otu-vs-asv>`_

**algorithms**

There are three general types of algorithm for clustering metabarcoding reads into OTU's or ASV's:

`Alignment based strategies <1\. Alignment-based strategy>`_
`De novo clustering - threshold <2\. De novo Clustering strategy - defined threshold_>`_
`De novo clustering - no threshold <3\. Clustering with guided clustering instead of thresholds>`_

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

Been in use for a long time so this process is well understood

There are more parameters than alignment strategy so the process is more complicated

3\. Clustering with guided clustering instead of thresholds
-----------------------------------------------------------

Not as arbitrary as threshold-based analysis
Accounts for sequencing errors

* `USEARCH <http://www.drive5.com/usearch/>`_
* `SWARM <https://github.com/torognes/swarm>`_
* `DADA2 <https://benjjneb.github.io/dada2/>`_

Identifying Chimeras
^^^^^^^^^^^^^^^^^^^^
[fill in]

[notes
1.  What proportion of the reads align to the reference?
2.  Chimera could be 2 species you haven't seen before
3.  Check OTU's individually
4.  Check against reference
5.  More abundant OTU's more likely to be real
6.  OTU's for every library and them split]


Taxonomic assignment:
^^^^^^^^^^^^^^^^^^^^^
[fill in]

Useful databases
----------------

Choosing a database ... [fill in]

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

Recommended packages and pipelines:
-----------------------------------

Depending on your skill level and needs for a given project the following packages and pipelines are Recommended:

# Table showing which processes each does and the environment needed e.g. galaxy. command line or R.

#
DADA2
Dadaist2
FROGS
PhyloSeq
Microbiome R package
MicrobiomeAnalyst
Rhea
Indecspecies
