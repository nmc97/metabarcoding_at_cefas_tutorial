Frogs
=====

`FROGS <https://github.com/geraldinepascal/FROGS#installation>`_

.. contents::
   :local:

Background & Installation
-------------------------

Frogs is a comprehensive pipeline for metabarcoding work. Users can run it within the command line or use galaxy to launch a graphical user interface.
It is designed to work regradless of the level of overlap of forward and reverse reads.

Installation on commandline:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Install using mamba:
This requires first downloading a file which specifies dependencies : `frogs-conda-requirements.yaml <https://github.com/geraldinepascal/FROGS/blob/master/frogs-conda-requirements.yaml>`_

The file contains something like this [check if it has been updated before using]:

.. code:: [yaml]

  channels:
    - conda-forge
    - bioconda
  dependencies:
  # bioconda
    - frogs =3.2.3
    - emboss =6.6
    - flash =1.2
    # need to be >=2.8
    - cutadapt =2.10
    # need to be >=2.1
    - swarm =3.0.0
    # need to be >= 2.17
    - vsearch =2.17.0
    - itsx =1.1.2
    - blast =2.10
    # - rdptool=2.0.3 # is already included in the frogs dependency
    - mafft =7.407
    - fasttree =2.1.9
    - bioconductor-phyloseq =1.34
    - bioconductor-deseq2 =1.30.1
  # conda-forge
    - r-base =4.0.5
    - r-essentials
    - r-phangorn =2.7.0
    - r-optparse
    - r-formattable
    - r-dt
    - r-plotly
    - r-gridextra

  ## packages included in previous package (possibly in multiple previous package)
  #   r-essentials includes:
      # pandoc
      # r-rmarkdown
      # r-dplyr
  #   bioconductor-phyloseq includes:
      # r-ggplot2
      # r-vegan
      # r-ape
      # r-scales
      # r-reshape2

Once you have saved the frogs-conda-requirements.yaml, run these commands:

.. code::

  mamba env create -n frogs@3.2.3 --file frogs-conda-requirements.yaml
  # activate your environment
  conda activate frogs@3.2.3

Test installation:

find the location of the frogs installaion and run test.sh

..code::







---
Author: Nicola Coyle
25/01/2022
