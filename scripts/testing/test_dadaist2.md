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


```
conda activate dadaist2 # not sure y the other one didn't work

cd /home/nc07/projects/metabarcoding/programs/dada2/dadaist2

# make a metadata file if one has not already been made
dadaist2-metadata -i /home/nc07/projects/metabarcoding/programs/dada2/dadaist2/data/16S_rename -o /home/nc07/projects/metabarcoding/programs/dada2/dadaist2/data/16S_rename/metadatafile.tsv

# main command - check parameters
dadaist2 \
-input-directory /home/nc07/projects/metabarcoding/programs/dada2/dadaist2/data/16S_rename \
-output-directory /home/nc07/projects/metabarcoding/programs/dada2/dadaist2/dadaist2-14-6_01 \
-database /home/nc07/projects/metabarcoding/programs/dada2/silva_nr99_v138.1_train_set.fa.gz \
-metadata /home/nc07/projects/metabarcoding/programs/dada2/dadaist2/data/16S_rename/metadatafile.tsv \
-threads 8 \
-trunc-len-1 250 \
-trunc-len-2 0 \
-min-qual 28 \
-maxee1 2 \
-maxee2 2 \
-save-rds \
-fastp \
-verbose

# export to get MicrobiomeAnalyst compatable files
dadaist2-exporter -i /home/nc07/projects/metabarcoding/programs/dada2/dadaist2/dadaist2-14-6_01
# make a multiqc report
dadaist2-mqc-report -i /home/nc07/projects/metabarcoding/programs/dada2/dadaist2/dadaist2-14-6_01  -o /home/nc07/projects/metabarcoding/programs/dada2/dadaist2/dadaist2-14-6_01/multiqc
# find alpha diversities - note -won't work if R packages aren't installed outside of a virtual env first Rhea automatically downloads before running if not available. Won't work in a virtual env with no internet. Missing GUniFrac for following script
dadaist2-normalize  -i /home/nc07/projects/metabarcoding/programs/dada2/dadaist2/dadaist2-14-6_01/feature-table.tsv -o /home/nc07/projects/metabarcoding/programs/dada2/dadaist2/dadaist2-14-6_01/normalise
dadaist2-alpha  -o /home/nc07/projects/metabarcoding/programs/dada2/dadaist2/dadaist2-14-6_01/alpha -i /home/nc07/projects/metabarcoding/programs/dada2/dadaist2/dadaist2-14-6_01/normalise/OTUs_Table-norm.tab # input file not found
```



test
/home/nc07/projects/gene_pools/data/GenePool_16S_MiSeq_DM13/16S_v3/

```
# User input: set file paths

in_dir=/home/user/path/to/read/directory/
out_dir=/home/user/path/to/read/directory/output
database=/home/nc07/path/to/database/silva_nr99_v138.1_train_set.fa.gz
meta=/home/user/path/to/metadatafile.csv # make one using dadaist2-metadata below

# example of test data
in_dir=/home/nc07/projects/gene_pools/data/GenePool_16S_MiSeq_DM13/16S_v3
out_dir=/home/nc07/projects/gene_pools/outputs_test_dadaist2
database=/home/nc07/projects/metabarcoding/programs/dada2/silva_nr99_v138.1_train_set.fa.gz
meta=/home/nc07/projects/gene_pools/data/GenePool_16S_MiSeq_DM13/16S_v3/metadata2.tsv # make one using dadaist2-metadata below

# activate conda environment
conda activate dadaist2-2

# optionally move into project directory
cd $in_dir

# make a metadata file if one has not already been made
dadaist2-metadata -i $in_dir -o $meta

# main command - check parameters
# note - if primers not supplied switch on fastp trimming. It will skip trimming if primers are not supplied and cutadapt trimming is selected.
time dadaist2 \
-input-directory $in_dir  \
-output-directory $out_dir \
-database  $database \
-metadata $meta \
-threads 12 \
-trunc-len-1 250 \
-trunc-len-2 0 \
-min-qual 28 \
-maxee1 2 \
-maxee2 2 \
-save-rds \
-fastp \
-verbose

# export to get MetagenomeAnalyist compatable files
dadaist2-exporter -i $out_dir
# make a multiqc report
dadaist2-mqc-report  -i $out_dir  -o $out_dir/multiqc
# find alpha diversities
dadaist2-normalize  -i $out_dir/feature-table.tsv -o $out_dir/normalise
