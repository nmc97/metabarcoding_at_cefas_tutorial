## Set up

### 1. Make a directory to run a the analysis
```
mkdir ~/metabarcoding_ws # make a base directory - if you are comfortable with linux you can put it somewhere else
mkdir ~/metabarcoding_ws/data # set directory to save data
mkdir ~/metabarcoding_ws/outputs # set directory to save outputs
mkdir ~/metabarcoding_ws/db # optionally put databases into a directory within workshop directory
```

### 2. [Download reference databases](https://metabarcoding-at-cefas-tutorial.readthedocs.io/en/latest/dada2.html#download-silva-datasets-curated-for-dada2)

Download Silva dataset for the tutorial (~ 1 mins 30 seconds, 131 Mb):
```
cd metabarcoding_ws/db # set a directory to store the data
time wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1 # grabs the file from the internet and downloads into the current directory
mv silva_nr99_v138.1_train_set.fa.gz?download=1 silva_nr99_v138.1_train_set.fa.gz # renames the file to remove "?download=1"
```

Download Silva *species* dataset for the tutorial (~ 30 seconds, 76 M):
```
cd metabarcoding_ws/db # set a directory to store the data
time wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1 # grabs the file from the internet and downloads into the current directory
mv silva_species_assignment_v138.1.fa.gz?download=1 silva_species_assignment_v138.1.fa.gz # renames the file to remove "?download=1"
```
### 3. Download test data

### 4. Set up conda environments

Set up Fastqc env:
```
mamba create -n fastqc fastqc MultiQC
```

Set up dada2 env:
```
mamba create -n dada2 r-essentials # setup a new environment and install r-essentials
conda activate dada2 # activate the new environment
mamba  install bioconductor-dada2 # install dada2
```

Set up dadaist2 env:
```
mamba create -n dadaist2
conda activate dadaist2
mamba install -c conda-forge -c bioconda dadaist2
mamba install bioconductor-dada2=1.20
mamba install -c conda-forge pyyaml # optional: needed to run dadaist2-mqc-report

```
