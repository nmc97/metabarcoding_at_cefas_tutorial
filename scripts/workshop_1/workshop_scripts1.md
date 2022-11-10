Metabarcoding workshop
29th July 2022
Weymouth

Author: Nicola Coyle
Based in part on the Dada2 tutorial: https://benjjneb.github.io/dada2/tutorial.html
Further resources can be found here: https://metabarcoding-at-cefas-tutorial.readthedocs.io/en/latest/index.html

For a script to apply to your own data see: https://github.com/nmc97/metabarcoding_at_cefas_tutorial/blob/main/scripts/run_dada2.R

## 1 Set up

### 1.1 Make a directory to run a the analysis
```
mkdir ~/metabarcoding_ws
mkdir ~/metabarcoding_ws/data # set directory to save data
mkdir ~/metabarcoding_ws/data/16S
mkdir ~/metabarcoding_ws/outputs # set directory to save outputs
mkdir ~/metabarcoding_ws/db # optionally put databases into a directory within workshop directory
```

### 1.2 [Download reference databases](https://metabarcoding-at-cefas-tutorial.readthedocs.io/en/latest/dada2.html#download-silva-datasets-curated-for-dada2)

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
### 1.3 Download test data

Jamie McMurtrie, Shayma Alathari, Dominique L. Chaput, David Bass, Camerson Ghambi, Joseph Nagoli, Jérôme Delamare-Deboutteville, Chadag Vishnumurthy Mohan, Joanne Cable, Ben Temperton, Charles R. Tyler, *Relationships between pond water and tilapia skin microbiomes in aquaculture ponds in Malawi*, Aquaculture, Volume 558, 2022, 738367, ISSN 0044-8486, https://doi.org/10.1016/j.aquaculture.2022.738367.

NCBI project code: PRJEB46984

Code for the project is available here: https://github.com/jamiemcm/Malawi_Tilapia_Microbiomes.

```bash
# make and activate a conda env for sra-tools
mamba create -n sra-tools
conda activate sra-tools

mkdir ~/metabarcoding_ws/data/16S
mkdir ~/metabarcoding_ws/data/16S/fastq

cd ~/metabarcoding_ws/data/16S

# download list of sample names and accession
wget https://raw.githubusercontent.com/nmc97/metabarcoding_at_cefas_tutorial/main/scripts/workshop_1/data/names_16S.txt

# prefetch data - uses file names_16S.txt to find the data (accession: q) download them (prefetch) and then rename the files (name: p)
while read p q; do
echo $p
time prefetch $q
done < names_16S.txt

# extract fastq files
while read p q; do
echo $p
ls $q/$q.sra
time fasterq-dump --split-files $q -o $q/$p
done < names_16S.txt

# move files and zip
mv ~/metabarcoding_ws/16S/*/*.fastq ~/metabarcoding_ws/16S/fastq
gzip ~/metabarcoding_ws/16S/fastq/*.fastq
```

Or - take it from my POD

``` bash
# this will take my fastq folder containing read files and place it in your 16S directory
rsync -ravz /home/nc07/metabarcoding_ws/data/16S/fastq/ ~/metabarcoding_ws/data/16S/
```
### 1.4 Set up conda environments

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

### 1.5 Set up RStudio with packages for vizualisation/ statistics (optionally)

``` R
[I may update this during the day]

```

## 2. Getting started on the day

### 2.1 Start a tmux session or screen
```
tmux -n meta
```
### 2.2 Start an interactive session for the day:
```
msub -I -q S30 -l procs=8,walltime=6:00:00 # change depending on what you think you need.
```

## 3. Quality control

## 3.1 Run [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [Multiqc](https://www.bing.com/search?q=multiqc&cvid=41c889025cf2474eae8615972425b160&aqs=edge..69i57j0l8j69i11004.1847j0j4&FORM=ANAB01&PC=U531)

```bash
#Activate conda environment:
	conda activate fastqc

#Set paths
	indir=~/metabarcoding_ws/data/16S/fastq
	outdir=~/metabarcoding_ws/outputs

#Make a directory to store fastqc files
	mkdir $outdir/fastqc_out

#Run fastqc:
	fastqc $indir/*.fastq.gz -t 8 -o $outdir/fastqc_out
# runs fastqc on all fastq.gz n current directory. -t threads to use

#Run multiqc on fastqc output folder.
#(It automatically detects fastqc outputs)
#change `--title` if you wish

	multiqc $outdir/fastqc_out/* -o $outdir/fastqc_out/multiqc_fastqc --title fastqc

#Note: force interative if there are lots of files (it will tell you if it wrote flat files instead)
	multiqc $outdir/fastqc_out/* -o $outdir/fastqc_out/multiqc_fastqc --interactive --title fastqc_interactive
```

## 4. Run [dada2 script](https://github.com/nmc97/metabarcoding_at_cefas_tutorial/blob/main/scripts/run_dada2.R) or [looped dada2 script](https://github.com/nmc97/metabarcoding_at_cefas_tutorial/blob/main/scripts/run_dada2.R)

## 4.1 Activate dada2 environment and start an R session

```bash
conda activate dada2
R
```

## 4.2 Set up R:
```R
#=========================#
# load libraries
#=========================#

library(dada2)

#=========================#
# setup file paths
#=========================#

# set outpath
outpath <- "~/metabarcoding_ws/outputs/dada2/"
# set path to data directory and list files
path <- "~/metabarcoding_ws/data/16S/fastq"

# if statement makes output folder if doesn't exist
if (file.exists(outpath)) {
 cat("The folder already exists")
} else {
 dir.create(outpath)
}

# check files exist in read file directory
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE)) ##change_me if the file names are different to this structure
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE)) ##change_me if the file names are different to this structure
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

## 4.3 Filter and trim:

### 4.3.1 Dada2 filter and trim

```R

# for completeness - check read quality in dada2
# make a pdf of the read quality output plot_read_quality_F.pdf
pdf(file = paste(outpath,"/plot_read_quality_F.pdf",sep=""), 7, 7)
plotQualityProfile(fnFs[1:12])
dev.off()
# make a pdf of the read quality output plot_read_quality_R.pdf
pdf(file = paste(outpath,"/plot_read_quality_R.pdf",sep=""), 7, 7)
plotQualityProfile(fnRs[1:12])
dev.off()

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(outpath, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) # forward reads
names(filtFs) <- sample.names
filtRs <- file.path(outpath, "filtered", paste0(sample.names, "_R_filt.fastq.gz")) # reverse reads
names(filtRs) <- sample.names

# Action required: Change trunclen=c(xx,xx) to match what you want to truncate your forward and reverse reads to. Similarly edit trimLeft=c(xx,xx), maxEE=x and truncQ=x. ##change_me
# https://rdrr.io/bioc/dada2/man/filterAndTrim.html
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(10,10), truncLen=c(200,160),
              maxN=0, maxEE=c(2,2), truncQ=c(2,2), rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=TRUE
head(out) # check how many reads have been lost in filtering step
# save out file in case environment shuts down
write.table(out,paste(outpath,"/out_save_F1.tsv",sep=""),sep="\t",row.names=T)

# make a pdf of the read quality output plot_read_quality_F.pdf
pdf(file = paste(outpath,"/plot_read_quality_filtered_F.pdf",sep=""), 7, 7)
plotQualityProfile(filtFs[1:12])
dev.off()
# make a pdf of the read quality output plot_read_quality_R.pdf
pdf(file = paste(outpath,"/plot_read_quality_filtered_R.pdf",sep=""), 7, 7)
plotQualityProfile(filtRs[1:12])
dev.off()

# Sometimes files lose all reads at this step. If that happens you need to exclude them from the rest the analysis.
# if statement checks if samples need to be removed and following statements remove them
removed_samples<-row.names(as.data.frame(out)[(which(as.data.frame(out)[,2] ==0)),])
#removed_samples<-c(removed_samples,"ESPCTD16S47_1.fastq.gz","ESPCTD16S80_1.fastq.gz","ESPCTD16S42_1.fastq.gz","ESPCTD16S36_1.fastq.gz")
removed_samples.names <- sapply(strsplit(basename(removed_samples), "_"), `[`, 1) # check this works for your data
if(length(removed_samples) != 0){
  filtFs<-filtFs[(!(names(filtFs) %in% removed_samples.names))] # remove from list of forward filtered read files
  filtRs<-filtRs[(!(names(filtRs) %in% removed_samples.names))] # remove from list of reverse filtered read files
  sample.names<-sample.names[!(sample.names) %in% removed_samples.names]
  out<-out[which(!(row.names(out)%in% removed_samples)),]
}

```

Stop here and check the filtering step before moving forward!
Repeat with new parameters if needed.

### 4.3.2 Run [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [Multiqc](https://www.bing.com/search?q=multiqc&cvid=41c889025cf2474eae8615972425b160&aqs=edge..69i57j0l8j69i11004.1847j0j4&FORM=ANAB01&PC=U531)

#### Open a new command-line session:

```bash
# start another interactive environement in a different session
msub -I -q S30 -l procs=8,walltime=0:30:00

#Activate conda environment:
	conda activate fastqc

#Set paths
	indir=~/metabarcoding_ws/outputs/dada2/filtered
	outdir=~/metabarcoding_ws/outputs/dada2/

#Make a directory to store fastqc files
	mkdir $outdir/filtered_fastqc_out

#Run fastqc:
	fastqc $indir/*.fastq.gz -t 8 -o $outdir/filtered_fastqc_out
# runs fastqc on all fastq.gz n current directory. -t threads to use

#Run multiqc on fastqc output folder.
#(It automatically detects fastqc outputs)
#change `--title` if you wish

	multiqc $outdir/filtered_fastqc_out/* -o $outdir/filtered_fastqc_out/multiqc_fastqc --title fastqc

#Note: force interative if there are lots of files (it will tell you if it wrote flat files instead)
	multiqc $outdir/filtered_fastqc_out/* -o $outdir/filtered_fastqc_out/multiqc_fastqc --interactive --title fastqc_filtered_interactive
```

If happy to proceed, calculate error rates:

## 4.4 Error rates
#### Back in R:
``` R
# calculate error rates for forward and reverse reads
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# make a pdf of the error output plot_errors.pdf
pdf(file = paste(outpath,"/plot_errors_F.pdf",sep=""), 7, 7)
plotErrors(errF, nominalQ=TRUE)
dev.off()
pdf(file = paste(outpath,"/plot_errors_R.pdf",sep=""), 7, 7)
plotErrors(errR, nominalQ=TRUE)
dev.off()
```

## 4.5 Dereplicate reads
``` R
derep_forward <- derepFastq(filtFs, verbose=TRUE)
# name the derep-class objects by the sample names
names(derep_forward) <- sample.names

derep_reverse <- derepFastq(filtRs, verbose=TRUE)
# name the derep-class objects by the sample names
names(derep_reverse) <- sample.names
```

## 4.6 Sample Inference
``` R
# Forward - this clusters forward reads into ASV's. later you will merge forward and reverse reads
dadaFs <- dada(derep_forward, err=errF, multithread=TRUE)
# Reverse
dadaRs <- dada(derep_reverse, err=errR, multithread=TRUE)

# merged paired read data
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# amplicon sequence variant table (ASV)
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # check size of table
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

```
## 4.7 remove chimeras
``` R
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) # check size of data.frame
print(paste(round(sum(seqtab.nochim)/sum(seqtab),4),"% of ASV's remain fter chimera removal: remaining ASV's:",dim(seqtab.nochim[2]),sep="")) # check proportion of reads are left after chimeras removed

# Save seqtab.nochim
saveRDS(seqtab.nochim, paste(outpath,"/seqtab_nochim.rds",sep=""))
```

## 4.8 Crude removal of negative controls
``` R
negative_controls<-c("Ts29")
seqtab.nochim.c<-seqtab.nochim
for(f in negative_controls){
  print(f)
  to_remove<-dim(seqtab.nochim.c[,which(seqtab.nochim.c[row.names(seqtab.nochim.c)==f,]>0)])[2]
  seqtab.nochim.c<-seqtab.nochim.c[,which(seqtab.nochim.c[row.names(seqtab.nochim.c)==f,]<=0)]
  print(paste(dim(seqtab.nochim.c)[2]," ASV's remain after removing", to_remove," ASV's present in negative controls:", f))
}
write.table(seqtab.nochim.c, file = paste(outpath,"/asv_table.dada2.tsv",sep=""),sep="\t",row.names=T)

# seqtab.nochim<-seqtab.nochim.c
# rm(seqtab.nochim.c)

```

## 4.9 track reads and ASV inference
``` R
# use this to determine which samples would be worth keeping in merged run

# set getN function
getN <- function(x) sum(getUniques(x))
# If you removed samples above you need to remove them again here:
# out<-out[which(row.names(out)!="sample_name1" & row.names(out)!="sample_name2" ),] ##change_me

# Combine filtering, dada2, and chimera removal data into one dataframe
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), rowSums(seqtab.nochim.c))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim","neg_removed")
rownames(track) <- sample.names
head(track) # check output

# write to file:
write.table(track, paste(outpath,"/track_reads.tsv", sep=""), sep="\t")
```

## 4.10 Taxonomy
``` R
# decide what reference dataset to use
# examples:

# dataset - SILVA:
taxa_with_bootstraps <- assignTaxonomy(seqtab.nochim, "~/metabarcoding_ws/db/silva_nr99_v138.1_train_set.fa.gz", outputBootstraps=TRUE, multithread=TRUE) ##change_me
taxa_all <- taxa_with_bootstraps$tax
# optionally add species:
taxa_species <- addSpecies(taxa_with_bootstraps$tax, "~/metabarcoding_ws/db/silva_species_assignment_v138.1.fa.gz")

# write to files
write.table(taxa_with_bootstraps, file = paste(outpath,"/taxa_with_bootstraps.tsv",sep=""),sep="\t",row.names=T)
write.table(taxa_all, file = paste(outpath,"/taxa.tsv",sep=""),sep="\t",row.names=T)
write.table(taxa_species, file = paste(outpath,"/taxa_species.tsv",sep=""),sep="\t",row.names=T)
```

Now you have completed the Dada2 aspect and have been left with output files:
ASV file:
Taxonomy files:
Track reads file:
Saved out file (can be discarded now):
Error rate plot pdfs:

## 5. Alternatively run [dadaist2]()

#### In the command-line:
```
msub -I -q S30 -l procs=8,walltime=1:00:00 # change depending on what you think you need.

in_dir=~/metabarcoding_ws/data/16S/fastq
out_dir=~/metabarcoding_ws/outputs/dadaist2
database=~/metabarcoding_ws/db/silva_nr99_v138.1_train_set.fa.gz
meta=~/metabarcoding_ws/outputs/dadaist2/metadatafile.csv # make one using dadaist2-metadata below
# activate conda environment
conda activate dadaist2

# optionally move into project directory
cd $in_dir
# make output directory
mkdir $out_dir

# make a metadata file if one has not already been made
dadaist2-metadata -i $in_dir -o $meta -1 _1 -2 _2

# main command - check parameters
# note - if primers not supplied switch on fastp trimming. It will skip trimming if primers are not supplied and cutadapt trimming is selected.
dadaist2 \
-input-directory $in_dir  \
-output-directory $out_dir \
-database  $database \
-metadata $meta \
-threads 8 \
-trunc-len-1 200 \
-trunc-len-2 160 \
-s1 0 \
-s2 0 \
-min-qual 28 \
-maxee1 2 \
-maxee2 2 \
-save-rds \
-1 _1 \
-2 _2 \
-verbose

# export to get MetagenomeAnalyist compatable files
dadaist2-exporter -i $out_dir
# make a multiqc report
dadaist2-mqc-report  -i $out_dir  -o $out_dir/multiqc
# find alpha diversities
dadaist2-normalize  -i $out_dir/feature-table.tsv -o $out_dir/normalise
# make a phyloseq object to import to R
dadaist2-phyloseqMake -i $out_dir
# use phyloseq object to automatically produce figures
dadaist2-taxaplot [options] -i phyloseq.rds -o $out_dir/plots
```

# 7. Remove contaminants - Decontam

Need a dataset with negative controls
https://github.com/benjjneb/decontam
https://bioconductor.org/packages/release/bioc/vignettes/decontam/inst/doc/decontam_intro.html
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("decontam")
```

Make a phyloseq object and then apply decontam functions.

# 6. Visualisation and statistics:

## 6.1 [MicrobiomeAnalyst](https://www.microbiomeanalyst.ca/)

To submit the files in the Dadaist MicrobiomeAnalyst directory you will first need to edit some of them.

The table.txt file needs to have file extensions removed from the sample names (use find and replace).

The taxa file doesn't contain enough delimiters. open in excel and save as csv file.

Download metadata file:

```
wget https://raw.githubusercontent.com/nmc97/metabarcoding_at_cefas_tutorial/main/scripts/workshop_1/data/metadata.csv
```

## 6.2 Phyloseq

For this you need R. It would be best to set up r studio on your computer and run this there.

You will need to create a Phyloseq object from the outputs of Dada2 or Dadaist2

For example: Dadaist two code for plots: https://quadram-institute-bioscience.github.io/dadaist2/notes/plot.html

## 6.3 Microbiome R package

For this you need R. It would be best to set up r studio on your computer and run this there.

You will need to create a Phyloseq object from the outputs of Dada2 or Dadaist2

https://microbiome.github.io/tutorials/

### 6.2.3 Installation

```
library(BiocManager)
BiocManager::install("microbiome")
```

Or

```
library(devtools) # Load the devtools package
install_github("microbiome/microbiome") # Install the package
```
