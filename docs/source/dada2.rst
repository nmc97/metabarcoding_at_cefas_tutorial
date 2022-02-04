Dada2
=====
`Dada2 - Divisive Amplicon Denoising Algorithm <https://github.com/benjjneb/dada2>`_

.. contents::
   :local:

Background & Installation
^^^^^^^^^^^^^^^^^^^^^^^^^

Dada2 is a widely used software for identifying ASV's in metabarcoding studies.

Papers:

`DADA2: High resolution sample inference from Illumina amplicon data, Challahan et al. 2016, Nature Methods <DADA2: High resolution sample inference from Illumina amplicon data>`_

`Bioconductor Workflow for Microbiome Data Analysis: from raw reads to community analyses, Challahan et al. 2016, F1000Research <https://f1000research.com/articles/5-1492`>_ (includes code)

`Exact sequence variants should replace operational taxonomic units in marker-gene data analysis, Challahan et al. 2017, Nature <https://www.nature.com/articles/ismej2017119>`_

`High-throughput amplicon sequencing of the full-length 16S rRNA gene with single-nucleotide resolution, Challahan et al. 2019, Nucleic acids rsearch 2019 <https://academic.oup.com/nar/article/47/18/e103/5527971>`_

Installation with mamba
-----------------------

Dada2 is an R package. One option for installing Dada2 in linux is to build a `Mamba <https://mamba.readthedocs.io/en/latest/user_guide/mamba.html>`_ environment. Some R packages can be install directly using `Mamba <https://mamba.readthedocs.io/en/latest/user_guide/mamba.html>`_. Others will need to be intstalled within R.

Create a Mamba environment:

.. code::

	mamba create -n dada2 r-essentials # setup a new environment and install r-essentials
	conda activate dada2 # activate the new environment
	mamba  install bioconductor-dada2 # install dada2

Alterantively - installing `Dada2` within R:

.. code::

	if (!requireNamespace("BiocManager", quietly = TRUE))
		install.packages("BiocManager")

	BiocManager::install("dada2")
	# note: not tested by NC

Dada2 Tutorial
^^^^^^^^^^^^^^

Dada2 tutorial can be found `here <https://benjjneb.github.io/dada2/tutorial.html>`_ : https://benjjneb.github.io/dada2/tutorial.html

Tutorial Data:
https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip

.. code::

	# download the tutorial data
	cd /path/to/tutorial_data/directory # set a directory to store the tutorial data
	wget https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip

Download Silva datasets. Curated for Dada2:
-------------------------------------------

Training taxonomy dataset for the tutorial. File location: https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1

.. code::

	# download silva dataset for the tutorial
	cd /path/to/data/directory # set a directory to store the data
	wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1 # grabs the file from the internet and downloads into the current directory
	mv silva_nr99_v138.1_train_set.fa.gz?download=1 silva_nr99_v138.1_train_set.fa.gz # renames the file to remove "?download=1"


Species Dataset. File location: https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1


.. code::

	# download silva species dataset for the tutorial
	cd /path/to/data/directory # set a directory to store the data
	wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1 # grabs the file from the internet and downloads into the current directory
	mv silva_species_assignment_v138.1.fa.gz?download=1 silva_species_assignment_v138.1.fa.gz # renames the file to remove "?download=1"


.. note:: Check that these downloadable datasets are up to date.

Follow up Statistics:
^^^^^^^^^^^^^^^^^^^^^

Install `phyloseq` and `Biostrings` in R
----------------------------------------

.. code::

	# install phyloseq within R - biocLite not working anymore, instead using BiocManager

	if (!require("BiocManager", quietly = TRUE))
		install.packages("BiocManager")

	BiocManager::install("phyloseq")

	# installing Biostrings

	if (!requireNamespace("BiocManager", quietly = TRUE))
		install.packages("BiocManager")

	BiocManager::install("Biostrings")

Dadaist2
^^^^^^^^

New : Wrapper - `Dadaist2: highway to R <https://quadram-institute-bioscience.github.io/dadaist2/>`_

Installation
------------

.. code ::

  mamba create -n dadaist2
  conda activate dadaist2
  mamba install -y -c conda-forge -c bioconda dadaist2

  # additionally install from github
  git clone https://github.com/quadram-institute-bioscience/dadaist2

Usage
-----

note 1 - file names must not start with a number

note 2 - can be run in POD using singularity and nextflow

Tutorial: https://quadram-institute-bioscience.github.io/dadaist2/tutorial

.. code ::

  dadaist2 -i data/16S/ -o example-output -d refs/SILVA_SSU_r138_2019.RData -t 8 -m metadata.tsv

From tutorial:

.. raw ::

  Briefly:

  -i points to the input directory containing paired end reads (by default recognised by _R1 and _R2 tags, but this can be customised)
  -o is the output directory
  -d is the reference database in DADA2 or DECIPHER format (we downloaded a DECIPHER database)
  -m link to the metadata file (if not supplied a blank one will be generated and used)
  -t is the number of processing threads

---
Author: Nicola Coyle
25/01/2022
