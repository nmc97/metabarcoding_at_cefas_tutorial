Resources
=========

+----------------------+--------------+-------------------------------------------+---------------------------------------------+-----------------------+
| Package/Pipeline     | How to use   | Processes                                 | features                                    | Issues                |
+----------------------+--------------+-------------------------------------------+---------------------------------------------+-----------------------+
| DADA2                | R            | Filtering and trimming reads              | Uses error model to cluster                 |                       |
|                      |              | Deduplicating reads                       |                                             |                       |
|                      |              | Clustering                                |                                             |                       |
|                      |              | Chimera removal                           |                                             |                       |
|                      |              | Taxonomic classification                  |                                             |                       |
+----------------------+--------------+-------------------------------------------+---------------------------------------------+-----------------------+
| Dadaist2             | Command-line | QC                                        | Run DADA2 using command-line with less code | Must use paired reads |
|                      |              | Filtering and trimming reads              |                                             |                       |
|                      |              |                                           |                                             |                       |
+----------------------+--------------+-------------------------------------------+---------------------------------------------+-----------------------+
| Phyloseq             | R            |                                           |                                             |                       |
+----------------------+--------------+-------------------------------------------+---------------------------------------------+-----------------------+
| MicrobiomeAnalyist   | Browser      | Visualise data, view statistical metrics  |                                             |                       |
+----------------------+--------------+-------------------------------------------+---------------------------------------------+-----------------------+
| Microbiome R package | R            | Visualise data, view statistical metrics  |                                             |                       |
+----------------------+--------------+-------------------------------------------+---------------------------------------------+-----------------------+
| Fastp                | Command-line | Read trimming and filtering               |                                             |                       |
+----------------------+--------------+-------------------------------------------+---------------------------------------------+-----------------------+
| Fastqc & multiqc     | Command-line | Quality assessment                        |                                             |                       |
+----------------------+--------------+-------------------------------------------+---------------------------------------------+-----------------------+
| Frogs                | Command-line |                                           |                                             |                       |
+----------------------+--------------+-------------------------------------------+---------------------------------------------+-----------------------+
| Frogs                | Galaxy       |                                           |                                             |                       |
+----------------------+--------------+-------------------------------------------+---------------------------------------------+-----------------------+


Visually explore data and statistics
------------------------------------

`MicrobiomeAnalyst <https://www.microbiomeanalyst.ca/>`_ :

Paper: `Using MicrobiomeAnalyst for comprehensive statistical, functional, and meta-analysis of microbiome data <https://www.nature.com/articles/s41596-019-0264-1>`_ `PDF <https://edisciplinas.usp.br/pluginfile.php/5269697/mod_resource/content/2/2020-Using%20MicrobiomeAnalyst%20for%20comprehensive%20statistical%2C%20functional%2C%20and%20meta-analysis%20of%20microbiome%20data.pdf>`_

Tutorial: `Metabarcoding 16s tutorial <https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/tutorials/MDP.pdf>`_

Uses: `Metacoder: An R package for visualization and manipulation of community taxonomic diversity data <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005404>`_ . Vizualising abundance comparisons in unique heatmap tree format. This is a good alternative to the stacked bar charts normally used to vizualise species abundance and compare conditions.

Where to find in MicrobiomeAnalyst: `Data Upload > Data Inspection > Data Filter > Normalization > Analysis Overview > Heat Tree`

.. note ::

  When to use:

  Excellent when you don't want to use R, and do not need to automate the output.

  Can be used directly after using Dadaist2 - outputs are supplied.

  Excellent for a first look at data.

  The drawbacks are that you cannot automate the process and you loose some control over the details.
  However, you can download the R scripts to ensure you can trace back what you have done.

Pipelines
---------

`micca <https://micca.readthedocs.io/en/latest/index.html>`_

"micca (MICrobial Community Analysis) is a software pipeline for the processing of amplicon sequencing data, from raw sequences to OTU tables, taxonomy classification and phylogenetic tree inference. The pipeline can be applied to a range of highly conserved genes/spacers, such as 16S rRNA gene, Internal Transcribed Spacer (ITS) 18S and 28S rRNA. micca is an open-source, GPLv3-licensed software."

`nf-core ampliseq <https://github.com/nf-core/ampliseq>`_

"nfcore/ampliseq is a bioinformatics analysis pipeline used for amplicon sequencing, supporting denoising of any amplicon and, currently, taxonomic assignment of 16S, ITS and 18S amplicons. Supported is paired-end Illumina or single-end Illumina, PacBio and IonTorrent data. Default is the analysis of 16S rRNA gene amplicons sequenced paired-end with Illumina."

Statistical analysis
--------------------

`PhyloSeq <https://micca.readthedocs.io/en/latest/phyloseq.html>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`<Shiny-phyloseq: Web application for interactive microbiome analysis with provenance tracking <https://academic.oup.com/bioinformatics/article/31/2/282/2365643>`_


`Rhea <https://lagkouvardos.github.io/Rhea/>`_ - statistical methods package for R
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Paper: `Rhea: a transparent and modular R pipeline for microbial profiling based on 16S rRNA gene amplicons <https://doi.org/10.7717/peerj.2836>`_

Installation :  Download a an extract folder into project directory - do this individually for each project. Use the installation script "install_packages.R" to install packages required (you can do this part once per R installation).


.. note ::

  When to use:

  If R is relatively new to you, this is designed to be simple to run and follow.


`Microbiome R packages <https://microbiome.github.io/tutorials/>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Paper:

Installation:

.. code ::

  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  # The following initializes usage of Bioc devel
  BiocManager::install(version='devel')

  BiocManager::install("microbiome")

Tutorial: `tutorial pages <https://microbiome.github.io/tutorials/>`_

.. note ::

  When to use:

  This package is extensive and can use data structured for use in phyloseq (a frequently used R package for diversity statistics).

  It appears to be in active use so bugs and issues will be addressed more easily by the community/ developers.

  The syntax is fairly clean and simple

  e.g.: to run an alpha diversity analysis looks like : `tab <-microbiome::alpha(pseq, index = "all")` vs  `tab <- estimate_richness(data) in phyloseq`)

  The microbiome R package will produce more alpha diversity meterics than pyloseq. This may be of use if you intend to use a different metric.

Exploring functional implications of community structures:
----------------------------------------------------------

`PICRUSt: Phylogenetic Investigation of Communities by Reconstruction of Unobserved States <http://picrust.github.io/picrust/>`_

Taxonomy
--------

`indicspecies: multivariate-analysis-indicator-value <https://www.rdocumentation.org/packages/indicspecies/versions/1.7.9/topics/indicspecies-package>`_

"This package provides a set of functions to assess the strength and statistical significance of the relationship between species occurrence/abundance and groups of sites. It is also possible to check the statistical significance of such associations."

Papers, chapters and commentary:
--------------------------------

Author = David Ryder & Nicola Coyle

Testing Alpha Diversity
^^^^^^^^^^^^^^^^^^^^^^^

`Comment in the Usearch Documentation <https://drive5.com/usearch/manual/alpha_diversity.html>`_

`Comment in the PhyloSeq FAQ <https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis>`_

`Paper discussing rarefaction of data <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531)>`_

`Paper discussion rarefaction of data <https://onlinelibrary.wiley.com/doi/epdf/10.1046/j.1461-0248.2001.00230.x>`_

`Chapter on species richness / alpha diversity metrics / population estimates 2001 <http://www.uvm.edu/~ngotelli/manuscriptpdfs/Chapter%204.pdf>`_

Testing Beta Diversity
^^^^^^^^^^^^^^^^^^^^^^

`Paper on normalisation prior to using beta diversity metrics <https://www.nature.com/articles/nmeth.2658>`_

Formats / standardisation
^^^^^^^^^^^^^^^^^^^^^^^^^

`Biom format <https://biom-format.org/documentation/biom_conversion.html>`_

Different algorithms
^^^^^^^^^^^^^^^^^^^^

`Dada2 Software <https://benjjneb.github.io/dada2/tutorial.html>`_

`Swarm Software <https://github.com/torognes/swarm>`_

`USearch Software <https://drive5.com/usearch/manual/uparse_pipeline.html>`_

Databases (lots of others)
^^^^^^^^^^^^^^^^^^^^^^^^^^

`PR2 database <https://github.com/pr2database/pr2database/releases>`_

`Silvia database <https://www.arb-silva.de/>`_

Fungi - ITS2
^^^^^^^^^^^^

`Best practices in metabarcoding of fungi: From experimental design to results <https://onlinelibrary.wiley.com/doi/10.1111/mec.16460#.YmMZICqe5zw.twitter>_`

`Database: Unite <https://unite.ut.ee/>`_

Nuclear ribosomal internal transcribed spacer (ITS)
Note: "By reanalysing published data sets, we demonstrate that operational taxonomic units (OTUs) outperform amplified sequence variants (ASVs) in recovering fungal diversity, a finding that is particularly evident for long markers. Additionally, analysis of the full-length ITS region allows more accurate taxonomic placement of fungi and other eukaryotes compared to the ITS2 subregion."

ANCOM
^^^^^

Analysis of composition of microbiomes: a novel method for studying microbial composition
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4450248/

---
Author: Nicola Coyle
25/01/2022
