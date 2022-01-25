Dada2
=====
`Divisive Amplicon Denoising Algorithm <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/>`_

.. contents::
   :local:
   
Background
^^^^^^^^^^
Dada2 is a widely used software for identifying ASV's in metabarcoding studies.

Installation
------------

Dada2 is an R package. One option for installing Dada2 in linux is to build a [Mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html) environment. Some R packages can be install directly using [Mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html). Others will need to be intstalled within R.

Create a Mamba env:
^^^^^^^^^^^^^^^^^^^

	mamba create -n dada2 r-essentials # setup a new environment and install r-essentials
	conda activate dada2 # activate the new environment
	mamba  install bioconductor-dada2 # install dada2


Install `phyloseq` and `Biostrings` in R
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code::

	# install phyloseq within R - biocLite not working anymore, instead using BiocManager

	if (!require("BiocManager", quietly = TRUE))
		install.packages("BiocManager")

	BiocManager::install("phyloseq")

	# installing Biostrings

	if (!requireNamespace("BiocManager", quietly = TRUE))
		install.packages("BiocManager")

	BiocManager::install("Biostrings")

Alterantively - installing `Dada2` within R:

.. code::

	# similarly you can install dada2 within R

	if (!requireNamespace("BiocManager", quietly = TRUE))
		install.packages("BiocManager")

	BiocManager::install("dada2")
	# not tested by NC

Dada2 Tutorial
--------------

Dada2 tutorial can be found `here <https://benjjneb.github.io/dada2/tutorial.html>`_ : https://benjjneb.github.io/dada2/tutorial.html

tutorial Data:
https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip

.. code::

	# download the tutorial data
	cd /path/to/tutorial_data/directory # set a directory to store the tutorial data
	wget https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip	

Download Silva datasets. Curated for Dada2:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

---
Author: Nicola Coyle
25/01/2022
