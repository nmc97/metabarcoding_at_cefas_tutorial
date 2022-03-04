Frogs
=====

`FROGS (Find Rapidly Otus Galaxy Solution) <https://github.com/geraldinepascal/FROGS#installation>`_

.. contents::
   :local:

Frogs - Background & Installation
---------------------------------

Frogs is a comprehensive pipeline for metabarcoding work. Users can run it within the command line or use galaxy to launch a graphical user interface.
It is designed to work regardless of the level of overlap of forward and reverse reads.

Installation on command-line:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Install using mamba:
This requires first downloading a file which specifies dependencies : `frogs-conda-requirements.yaml <https://github.com/geraldinepascal/FROGS/blob/master/frogs-conda-requirements.yaml>`_

The file contains something like this [check if it has been updated before using]:

.. code::

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

note: the publisher recommends making an environment called `frogs@3.2.3`. This will make file locations more difficult to access due to the `@` symbol. Best to use frogs_3.2.3 instead.

.. code::

  mamba env create -n frogs_3.2.3 --file frogs-conda-requirements.yaml
  # activate your environment
  conda activate frogs_3.2.3

Test installation:

find the location of the frogs installaion ( likely will look like:/home/$Username/mambaforge/envs/frogs_3.2.3/share/FROGS-3.2.3/) and run test.sh (within test folder)

.. code ::

  cd $path_to_frogs/test # for me path_to_frogs was in /home/nc07/mambaforge/envs/frogs@3.2.3/share/FROGS-3.2.3/test
  cd /home/nc07/mambaforge/envs/frogs@3.2.3/share/FROGS-3.2.3/test

  sh test.sh ../ 1 2 res

This resulted in an error because cutadapt has not been installed properly
Similarly, swarm and emboss needed to be installed individually. In addition openssl need to be downgraded.

.. code ::

  mamba install cutadapt=2.10
  mamba install swarm=3.10
  mamba install emboss=6.6.0
  mamba install openssl=1.1.1m

If you encounter a different error, try to run the command that the program failed on. take a look at the output and see if you can decipher what went wrong. You can check the installation of each program by typing it into the command line separately, if you know it's name.

Tutorial
--------

`tutorial <https://tutorials.migale.inra.fr/posts/frogs-16s/#:~:text=FROGS%20%5B%201%5D%20is%20a%20tool%20dedicated%20to,performed%20on%20the%20Migale%20cluster%20migale.jouy.inrae.fr%20and%20rstudio.migale.inrae.fr.>`_

Parameter Considerations
------------------------

SWARM - d
^^^^^^^^^

SWARM uses a parameter d to determine which sequences to add to a cluster. If a new sequence is d distance from an sequences in an existing cluster, the new sequence is added to this cluser. `d` is 1 by default.

In `Ershova et al. 2021 <https://academic.oup.com/icesjms/article/78/9/3342/6360557#323435484>`_ they use `d = 13`:

"Step-by-step clustering was performed in SWARM 2.1.13 (Mahé et al., 2015) using a distance value of d = 13 to cluster individual sequences into molecular operational taxonomic units (MOTUs). This distance value has previously been used to cluster similar datasets using the same COI fragment (e.g. Bakker et al., 2019; Antich et al., 2020; Atienza et al., 2020)."

---
Author: Nicola Coyle
25/01/2022
