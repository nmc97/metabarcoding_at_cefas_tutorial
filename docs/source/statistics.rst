==========
Statistics
==========

Abundance statistics of metabarcoding data

Normalisation
^^^^^^^^^^^^^

Do not remove singletons before calculating chao1 predicted diversity. This process relies on knowing the number of singleton and doubleton OTU's/ASV's. This also means that is there are many spurious singleton/ doubleton OTU's (due to sequencing error etc.), then the chao estimates will be poor.

To rarefy or to not rarefy
^^^^^^^^^^^^^^^^^^^^^^^^^^

Rarefying data is a method of normalising datasets so that they can be directly compared without biases.
The main bias being avoided in metabarcoding is read depth -
where some samples have many more sequences than others the diversity they will uncover will surely bee different purely because of these extra reads.

The procedure followed is normally

  - choose a minimum number of reads that must be present in a sample

  - remove samples with fewer reads than this chosen threshold

  - randomly Subsample all samples until they have the same number of reads

In phyloseq use `rarefy_even_depth`.

The major drawback of this approach is loss of data. "Despite its current popularity in microbiome analyses rarefying biological count data is statistically inadmissible because it requires the omission of available valid data." `(ref) <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531>`_

Additionally,

`Benchmarking microbiome transformations favors experimental quantitative approaches to address compositionality and sampling depth biases <https://www.nature.com/articles/s41467-021-23821-6>`_

Alpha Diversity
^^^^^^^^^^^^^^^

:term:`Alpha diversity`


Beta Diversity
^^^^^^^^^^^^^^

:term:`Beta diversity`

:term:`Bray Curtis distance`

Differential Abundance
^^^^^^^^^^^^^^^^^^^^^^

`Analysis of microbial compositions: a review of normalization and differential abundance analysis <https://www.nature.com/articles/s41522-020-00160-w>`_

---
Author: Nicola Coyle
25/01/2022
