# Test dadaist2

https://quadram-institute-bioscience.github.io/dadaist2/tutorial

```

cd /home/nc07/projects/metabarcoding/programs/dada2/dadaist2/

conda activate dadaist2

dadaist2-metadata -i data/16S > metadata.tsv

dadaist2 -i data/16S/ -o example-output -d refs/SILVA_SSU_r138_2019.RData -t 8 -m metadata.tsv

```



**Error: **

`Error in validObject(.Object) :
  invalid class “SRFilterResult” object: superclass "Mnumeric" not defined in the environment of the object's class`

dada2 `version 1.18.0`

**similar error:** https://github.com/benjjneb/dada2/issues/1378?msclkid=22e9afc7c56511ecb224dee2d9d45399


## Reinstall dadaist2 and upgrade dada2

Set up new environment and test dada2 v 1.20 with dadaist2

``` bash
mamba create -n dadaist2-2
conda activate dadaist2-2
mamba install -y -c conda-forge -c bioconda dadaist2 # takes a while!!
mamba install bioconductor-dada2=1.20
```

**additionally install from github**

``` bash
cd /path/to/conda/environment/dadaist2/directory
git clone https://github.com/quadram-institute-bioscience/dadaist2
```

# why did I do this last tie ??not replacted
replace r curl version and reinstall PhyloSeq
Rcurl Version: 1.98-1.5

``` R
uninstall.packages(RCurl)
```

## Test:

``` bash

cd
conda activate dadaist2-2
dadaist2 -i data/16S/ -o example-output -d refs/SILVA_SSU_r138_2019.RData -t 8 -m metadata.tsv

dadaist2-exporter -i example-output


```
