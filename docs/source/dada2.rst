Dada2
=====
`Dada2 - Divisive Amplicon Denoising Algorithm <https://github.com/benjjneb/dada2>`_

.. contents::
   :local:

Dada2 - Background & Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Dada2 is a widely used software for identifying ASV's in metabarcoding studies.

Papers:

`DADA2: High resolution sample inference from Illumina amplicon data, Challahan et al. 2016, Nature Methods <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/>`_

`Bioconductor Workflow for Microbiome Data Analysis: from raw reads to community analyses, Challahan et al. 2016, F1000Research <https://f1000research.com/articles/5-1492>`_ (includes code)

`Exact sequence variants should replace operational taxonomic units in marker-gene data analysis, Challahan et al. 2017, Nature <https://www.nature.com/articles/ismej2017119>`_

`High-throughput amplicon sequencing of the full-length 16S rRNA gene with single-nucleotide resolution, Challahan et al. 2019, Nucleic acids rsearch 2019 <https://academic.oup.com/nar/article/47/18/e103/5527971>`_

Installation with mamba
-----------------------

Dada2 is an R package. One option for installing Dada2 in Linux is to build a `Mamba <https://mamba.readthedocs.io/en/latest/user_guide/mamba.html>`_ environment. Some R packages can be install directly using `Mamba <https://mamba.readthedocs.io/en/latest/user_guide/mamba.html>`_. Others will need to be installed within R.

Create a Mamba environment:

.. code::

	mamba create -n dada2 r-essentials # setup a new environment and install r-essentials
	conda activate dada2 # activate the new environment
	mamba  install bioconductor-dada2 # install dada2

Installation within R:
--------------------------------------------

.. code::

	if (!requireNamespace("BiocManager", quietly = TRUE))
		install.packages("BiocManager")

	BiocManager::install("dada2")


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


.. note:: Check that these downloadable datasets are up to date before downloading!

DADA2 Usage
^^^^^^^^^^^

**Optional:** Before using Dada2 you can do some quality testing suing fastqc and multiqc. See script `<https://github.com/nmc97/metabarcoding_at_cefas_tutorial/blob/main/scripts/metabarcoding_qc.sh>`_

To start using Dada2 on POD, start an interactive session. Navigate to a working directory, activate dada2 conda environment, and start R.

.. code::

  msub -I -q S30 -l procs=12,walltime=6:00:00 # change depending on what yo think you need.

  conda activate dada2

  cd /path/to/project/folder

  # start R
  R

Once in R you can follow the script `run_dada2.R <https://github.com/nmc97/metabarcoding_at_cefas_tutorial/blob/main/scripts/run_dada2.R>`_.

Dada2 Tutorial
^^^^^^^^^^^^^^

To familiarise yourself with Dada2 see the Dada2 tutorial `here <https://benjjneb.github.io/dada2/tutorial.html>`_ : https://benjjneb.github.io/dada2/tutorial.html

The tutorial data is available here:
https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip

.. code::

	# download the tutorial data
	cd /path/to/tutorial_data/directory # set a directory to store the tutorial data
	wget https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip


Alternative tutorial:
^^^^^^^^^^^^^^^^^^^^^
https://replikation.github.io/bioinformatics_side/metagenome/metabarcoding/

Dadaist2
^^^^^^^^

Dadaist2 is a command line wrapper for Dada2

`Dadaist2: highway to R <https://quadram-institute-bioscience.github.io/dadaist2/>`_

.. note::

  When to use:

  If you like working within the command line instead of R, this could be ideal. Familiarity with Dada2 methods is necessary to ensure the parameters involved are correct for your data. It has many automatically generted outputs that may be very useful e.g. MicrobiomeAnalyst, phyloseq and Rhea input files, and very nice html log files. Rhea is incorporated into Dadaist2 so some statistical analysis can be conducted using this package.

Installation
------------

.. code ::

  mamba create -n dadaist2
  conda activate dadaist2
  mamba install -y -c conda-forge -c bioconda dadaist2
  mamba install bioconductor-dada2=1.20
  mamba install -c conda-forge pyyaml # optional: needed to run dadaist2-mqc-report

**Additionally install from github**
Navigate to the directory that has been made for the new environment:
`cd /path/to/conda/environment/dadaist2/directory`
eg:

.. code::

  $ whereis dadaist2 # finds where the command is located
  dadaist2: /home/username/mambaforge/envs/dadaist2/bin/dadaist2
  $ cd /home/username/mambaforge/envs/dadaist2/bin/ # navigate to that directory
  $ git clone https://github.com/quadram-institute-bioscience/dadaist2 # install from github

Install Rhea packages for downstream analysis. Rhea is used in some dadaist2 scripts to assess diversity. In order to use these scripts within a POD virtual environment which cannot access the internet to download new packages, you will need to down;oad Rhea prerequisites yourself first.

Open R and use the following to check if GUniFrac and vegan are installed and install them.

.. code::

  # code from https://github.com/Lagkouvardos/Rhea/blob/master/install_packages.R
  # Check if required packages are already installed, and install if missing
  packages <- c("GUniFrac","vegan")

  # Function to check whether the package is installed
  InsPack <- function(pack)
  {
    if ((pack %in% installed.packages()) == FALSE) {
      install.packages(pack,repos = "http://cloud.r-project.org/")
    }
  }

  # Applying the installation on the list of packages
  lapply(packages, InsPack)

  # Make the libraries
  lib <- lapply(packages, require, character.only = TRUE)

  # Check if it was possible to install all required libraries
  flag <- all(as.logical(lib))

Usage
-----

Note - File names must not start with a number! An unfortunate issue, but likely due to R not liking names beginning with a number.

You can follow a tutorial and view documentation here: https://quadram-institute-bioscience.github.io/dadaist2/tutorial. Note that the test data results cannot be loaded into MicrobiomeAnlaysist becasue there are too many OTU's unique to each sample, meaning they have nothing to show.
Download github code to access test data:

.. code ::

  git clone https://github.com/quadram-institute-bioscience/dadaist2
  cd dadaist2

Minimal use case:

.. code ::

  dadaist2 -i data/16S/ -o example-output -d refs/SILVA_SSU_r138_2019.RData -t 8 -m metadata.tsv

  # Briefly:

  # -i points to the input directory containing paired end reads (by default recognised by _R1 and _R2 tags, but this can be customised)
  # -o is the output directory
  # -d is the reference database in DADA2 or DECIPHER format (we downloaded a DECIPHER database)
  # -m link to the metadata file (if not supplied a blank one will be generated and used)
  # -t is the number of processing threads

More extensive example:

.. code::

  conda activate dadaist2 # not sure y the other one didn't work

  cd /home/user/path/to/project/directory/

  # make a metadata file if one has not already been made
  dadaist2-metadata -i /home/user/path/to/project/directory/ -o  /home/user/path/to/project/directory/metadatafile.tsv

  # main command - check parameters
  dadaist2 \
  -input-directory /home/user/path/to/read/directory/ \
  -output-directory /home/user/path/to/read/directory/output \
  -database /home/user/path/to/database/silva_nr99_v138.1_train_set.fa.gz \
  -metadata /home/user/path/to/metadatafile.csv \
  -threads 12 \
  -trunc-len-1 250 \
  -trunc-len-2 0 \
  -min-qual 28 \
  -maxee1 2 \
  -maxee2 2 \
  -save-rds \
  -verbose

  # export to get MetagenomeAnalyist compatable files
  dadaist2-exporter -i /home/user/path/to/read/directory/output
  # make a multiqc report
  dadaist2-mqc-report  -i /home/user/path/to/read/directory/output  -o /home/user/path/to/read/directory/output/multiqc
  # find alpha diversities
  dadaist2-normalize  -i /home/user/path/to/read/directory/output/MetagenomeAnalyist -o OUTDIR

You can follow the script `run_dadaist2.sh <https://github.com/nmc97/metabarcoding_at_cefas_tutorial/blob/main/scripts/run_dadaist2.sh>`_ to apply the above to your data with more ease.

Note : if primers not supplied switch on fastp trimming using the `--fastp` flag. It will skip trimming entirely if primer sequences are not supplied and the default cutadapt trimming is selected.


Plotting Taxonomy Dadaist2 vs PhyloSeq
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use script `dadaist2-taxplot` in Dadaist2

`Notes on comparison <https://quadram-institute-bioscience.github.io/dadaist2/notes/6_Rscripts.html>`_
`Phyloseq script <https://quadram-institute-bioscience.github.io/dadaist2/notes/plot.html>`_

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


---
Author: Nicola Coyle
25/01/2022
