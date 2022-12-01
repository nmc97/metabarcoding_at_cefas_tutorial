# Metabarcoding workshop
*29th July 2022
Weymouth*

**Author:** Nicola Coyle
Based in part on the Dada2 tutorial: https://benjjneb.github.io/dada2/tutorial.html

To learn more about Dada2 please see the [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/) and [manual](https://www.bioconductor.org/packages/devel/bioc/manuals/dada2/man/dada2.pdf)

Further resources can be found here: https://metabarcoding-at-cefas-tutorial.readthedocs.io/en/latest/index.html

If you would like to try this on your own data you can use the script supplied here: https://github.com/nmc97/metabarcoding_at_cefas_tutorial/blob/main/scripts/run_dada2.R

## 1. Set up

First, let's set up a directory to do our work in. We will call is metabarcoding_ws and place it in our home directory (`~/` or `/home/$USER`).

Next we will download reference databases for classifying our barcodes. in this example we will be using the SILVA database.

Then we will download our data and get started!

**Some of the basic commands we will use:**

`mkdir /path/to/directory`: This makes a directory

`cd /path/to/directory`: This changes directory

`mv /path/to/file /path/to/new_file_location`: This will move a file or folder from one place to another. You can use this to change the file name.

### 1.1 Make a directory to run a the analysis
``` bash
mkdir ~/metabarcoding_ws # makes a directory metabarcoding_ws in the home directory `~/`
mkdir ~/metabarcoding_ws/data # make a directory to save data
mkdir ~/metabarcoding_ws/data/16S # specify a data directory for 16S data
mkdir ~/metabarcoding_ws/outputs # set directory to save outputs
mkdir ~/metabarcoding_ws/db # optionally put databases into a directory within workshop directory
```

### 1.2 [Download reference databases](https://metabarcoding-at-cefas-tutorial.readthedocs.io/en/latest/dada2.html#download-silva-datasets-curated-for-dada2)

Download Silva dataset for the tutorial (~ 1 mins 30 seconds, 131 Mb):
```
# move to the directory where we will store the databases using `cd` (change directory)
cd metabarcoding_ws/db

# wget grabs the file from the internet and downloads into the current directory
wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1

# rename the file to remove "?download=1". Uses command `mv`
mv silva_nr99_v138.1_train_set.fa.gz?download=1 silva_nr99_v138.1_train_set.fa.gz
```

Download Silva *species* dataset for the tutorial (~ 30 seconds, 76 M):
```
cd metabarcoding_ws/db # set a directory to store the data
wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1 # grabs the file from the internet and downloads into the current directory
mv silva_species_assignment_v138.1.fa.gz?download=1 silva_species_assignment_v138.1.fa.gz # renames the file to remove "?download=1"
```

Check if this is the latest version of the dataset here: https://www.arb-silva.de/download/archive/

### 1.3 Download test data

Here we will use a dataset produced in collaboration with Cefas:

Jamie McMurtrie, Shayma Alathari, Dominique L. Chaput, David Bass, Camerson Ghambi, Joseph Nagoli, Jérôme Delamare-Deboutteville, Chadag Vishnumurthy Mohan, Joanne Cable, Ben Temperton, Charles R. Tyler, *Relationships between pond water and tilapia skin microbiomes in aquaculture ponds in Malawi*, Aquaculture, Volume 558, 2022, 738367, ISSN 0044-8486, https://doi.org/10.1016/j.aquaculture.2022.738367.

NCBI project code: PRJEB46984

Code for the project is available here: https://github.com/jamiemcm/Malawi_Tilapia_Microbiomes.

There are two ways to get the data added to your POD:

**Option 1: Take it from my POD (recommended)**

``` bash
# this will take my fastq folder containing read files and place it in your 16S directory
rsync -ravz /home/nc07/metabarcoding_ws/data/16S/fastq ~/metabarcoding_ws/data/16S/
```

**Option 2: Download from SRA (long)**
Uses [SRA-tools](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump)

This code requires a list of sample names and accession numbers to run. I have compile one here: https://raw.githubusercontent.com/nmc97/metabarcoding_at_cefas_tutorial/main/scripts/workshop_1/data/names_16S.txt

```bash
# make and activate a conda env for sra-tools
mamba create -n sra-tools
conda activate sra-tools

# Create directories to download files into
mkdir ~/metabarcoding_ws/data/16S
mkdir ~/metabarcoding_ws/data/16S/fastq

cd ~/metabarcoding_ws/data/16S

# download list of sample names and accession numbers
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

# move files from individual folders into ~/metabarcoding_ws/16S/fastq and gzip them
mv ~/metabarcoding_ws/16S/*/*.fastq ~/metabarcoding_ws/16S/fastq
gzip ~/metabarcoding_ws/16S/fastq/*.fastq
```

### 1.4 Set up conda environments

**Set up Fastqc env:**
``` bash
# create an environment called fastqc and install fastqc and multiqc
mamba create -n fastqc fastqc MultiQC
```

If this is your first time creating a conda environment, or it has been a while, here's a breakdown of what just happened:

"mamba": calls the program mamba
"create": tells mamba you want to create a new environment
"-n fastqc": names the new environment fastqc
"fastqc": installs the program named fastqc
"MultiQC": installs the program MultiQC

You can list as many programs as you want after naming the environment and Mamba will attempt to install them.

However, it is good practice to keep individual environments for specific programs. This is because when you add a new program to environment, it may require different dependencies than programs already in the environment. Thus, the dependencies could clash and cause one or many programs so stop working. Here we will keep the number of programs in each environment low, and try to group them by task. Fastqc and MultiQC are frequently used together and work well in and envrionment.

In the next environment we set up R using the package r-essentials. Then we can install dada2, which runs in R:

**Set up dada2 env:**
``` bash
mamba create -n dada2 r-essentials # setup a new environment and install r-essentials
conda activate dada2 # activate the new environment
mamba  install bioconductor-dada2 # install dada2
```

Lastly we will install dadasit2, which is a wrapper program for dada2. This means that they bundled up code for running dada2 into a scripts so that you can run dada2 using a single line of code in the command-line. This could be useful if you have data that you want to quickly analyse which isn't too complex.

**Set up dadaist2 env:**

``` bash
mamba create -n dadaist2
conda activate dadaist2
mamba install -c conda-forge -c bioconda dadaist2
mamba install bioconductor-dada2=1.20
mamba install -c conda-forge pyyaml # optional: needed to run dadaist2-mqc-report
```

## 2. Getting our POD session set up

### 2.1 Start a screen session

When you are logged in to POD, and your connection drops, you will loose all the progress you had made in that terminal.

Screen allows you to set up a session or a "screen" in which to do your work. You can shut the screen down and it will continue to run in the background. We will set one up for this session:

``` bash
screen -S meta
```

You can leave the screen by pressing `ctrl + a then d`

Restart the screen by typing `screen -r meta`

### 2.2 Start an interactive session for the day (using POD at Cefas):

Where we are all working now, there isn't much computing power available. In POD we will need to each ask for more by starting an interactive session.

Here, you can change the S30 to B30 if there is no more space left on S30.

See [here](https://cefas.sharepoint.com/sites/PODHighPerformanceComputingUsers/SitePages/Running%20Jobs%20on%20MT2.aspx) for more information.

``` bash
msub -I -q S30 -l procs=8,walltime=7:00:00 # change depending on what you think you need.
```

## 3. Quality control

First thing we need to do is check the quality of the data. We will run fastqc which will assess all the fastq files filled with reads. It will tell us how many reads we have and what they look like in regards to length and quality.

## 3.1 Run [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [Multiqc](https://www.bing.com/search?q=multiqc&cvid=41c889025cf2474eae8615972425b160&aqs=edge..69i57j0l8j69i11004.1847j0j4&FORM=ANAB01&PC=U531)

The first thing we will do is activate the conda envrionment we made earlier called "fastqc". This environment contains programs fastqc and MultiQC.

Then we will set some paths. This means that we will set two variables `$indir` and `$outdir` which contain the paths on the files system. Variables in bash are denoted using a `$`. The first path, `$indir`, will contain the path to the input data for the experiment. We put this data in he directory (folder) called `~/metabarcoding_ws/data/16S/fastq`, so this will be set as our, `$indir`, variable.

If the output directory for fastqc doesn't exist it will not be happy, so we need to make it. To make a directory we use the command `mkdir`.

After we make this directory we run fastqc.

Lastly we will run MultiQC, which collates all the outputs from fastqc and puts them into one report.

Okay, let's try it:

```bash
# Activate conda environment:
	conda activate fastqc

# Set paths
	indir=~/metabarcoding_ws/data/16S/fastq # this is where our input data lives
	outdir=~/metabarcoding_ws/outputs # this is where all our outputs will go

# Make a directory to store fastqc files
	mkdir $outdir/fastqc_out
```

Lets check that that worked by listing all the files in our output directory:

`ls $outdir`

`ls` is a command to list files. Here it will list all the files in `$outdir`.

If it gives `fastqc_ws` as an output on the terminal, then we successfully created that directory.

Now lets run FastQC:

``` bash
# Run fastqc:
	fastqc $indir/*.fastq.gz -t 8 -o $outdir/fastqc_out
# runs fastqc on all fastq.gz n current directory. -t threads to use

# Run multiqc on fastqc output folder $outdir/fastqc_out/ usng wild card * to catch all the fies in the directory
# (It automatically detects fastqc outputs)
# change `--title` if you wish

	multiqc $outdir/fastqc_out/* -o $outdir/fastqc_out/multiqc_fastqc --title fastqc
```

This gives a message at the end:
``` bash
|           multiqc | 8 flat-image plots used in the report due to large sample numbers
|           multiqc | To force interactive plots, use the '--interactive' flag. See the documentation.
```

So lets add the `--interactive` flag to the command, and change the output name:
``` bash
# force interative if there are lots of files (it will tell you if it wrote flat files instead)
	multiqc $outdir/fastqc_out/* -o $outdir/fastqc_out/multiqc_fastqc --interactive --title fastqc_interactive
```

Now we will check the quality of the data and decide how we would like to trim and filter them in the next step.

You can do this in MobaXterm by navigating to the MultiQC output directory and selecting the `.html` in the file.

If we would like to, we can find out the full path to this directory by using `echo`:

`echo $outdir/fastqc_out/multiqc_fastqc`: this will print the output of `$outdir/fastqc_out/multiqc_fastqc` to the command line with the variable `$outdir` filled in.

Once you find the `html` file, right click and select `open with`. Then select a browser.

Question: What do you notice about the reads?  

## 4. Run [dada2 script](https://github.com/nmc97/metabarcoding_at_cefas_tutorial/blob/main/scripts/run_dada2.R) or [looped dada2 script](https://github.com/nmc97/metabarcoding_at_cefas_tutorial/blob/main/scripts/run_dada2.R)

## 4.1 Activate dada2 environment and start an R session

Firstly, lets deactivate the fastqc environment now that we are done with it:

```bash
# deactivate fastqc conda environment
conda deactivate fastqc
```

Next we will activate out dada2 environment, and start an R session. We do this in the terminal by typing R.

```bash
# activate dada2 conda environment
conda activate dada2

# start R session
R
```

## 4.2 Set up R:

First thing to do in R is load the dada2 library `library(dada2)`. You will need to do this any time that you restart the R session. If you leave the session and return, loading the R session again, it will not reload the libraries for you.



```R
#=========================#
# load libraries
# we only need dada2
#=========================#

library(dada2)

```

Next we will do the same as above and set some variables (outpath, path). Note that in R we don't need to indicate something is a variable using `$` as in bash.

``` R
#=========================#
# setup file paths
#=========================#

# set path to output directory outpath
outpath <- "~/metabarcoding_ws/outputs/dada2/"
# set path to data directory
path <- "~/metabarcoding_ws/data/16S/fastq"

# check if the output folder exists.
# Use this if statement to make an output folder if doesn't exist:
if (file.exists(outpath)) {
 cat("The folder already exists")
} else {
 dir.create(outpath)
 cat(paste("making a directory:", outpath, "\n"))
}

# check input files exist in read file directory
# list.files will list all the files in path
list.files(path)
```

dada2 uses lists of file names, so we need to produce them. We will produce a list of forward read files `fnFs` and reverse reads `fnRs`.

> *Note:* if your files end in something different to `_1.fastq.gz` or `_2.fastq.gz` you will need to change the code to match your data.
> In particular: `pattern="_1.fastq"` and `pattern="_2.fastq"`
>
> For example if you're files have an extension `_R1_001.fastq.gz` and `R2_001.fastq.gz`, you could use the pattern=`_R1_001.fastq` and pattern=`_R2_001.fastq`.
> You can use any pattern that allows the function to find every forward read file and each reverse read file correctly 

``` R
# assuming Forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_R2_001.fastq

# list forward read files

fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE)) # this will find all files otaining "_1.fastq"

# list reverse read files

fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE)) # this will find all files otaining "_2.fastq"

# Extract and keep a list of sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

This last line may need your attention when it comes to your own data.

Firstly, make sure that you know what sample names you would like to apply to your data and that they appear in your file names. 

This line will extract the string that comes before the first underscore. If you would like it to extract everything before the second underscore you can change the number `1` to a `2`. 

For example:

`Sample_name_1.fastq.gz`

`sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)`

However, if you have samples that have underscores in the sample names, this line will not work properly. In this case it will only use the word before the first underscore "_".
For example, if your samples have “_” within the name, try to find a different point at which to split the filename up. 

For example:

`SAMPLE_NAME_16S_R1.fastq.gz`

`sample.names <- sapply(strsplit(basename(fnFs), "_16S"), `[`, 1)`

## 4.3 Filter and trim:

Now that we have set up all the file paths, we can have a second look at the quality of the reads. This time, using dada2.

Then we can decide how we would like to trim and filter them.

### 4.3.1 Dada2 filter and trim

For completeness - we will check read quality in dada2. dada2 has a function `plotQualityProfile` which can be used to view the read quality of your samples.

```R
# make a pdf of the read quality output plot_read_quality_F.pdf
# 12 samples included. The first 12 samples are included by using [1:12]
pdf(file = paste(outpath,"/plot_read_quality_F.pdf",sep=""), 7, 7)
plotQualityProfile(fnFs[1:12])
dev.off()

# make a pdf of the read quality output plot_read_quality_R.pdf
# 12 samples included

pdf(file = paste(outpath,"/plot_read_quality_R.pdf",sep=""), 7, 7)
plotQualityProfile(fnRs[1:12])
dev.off()
```

> Note: you will see the following error message:
`guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> = "none")` instead.
Ignore this.

>See the [dada2 manual](https://www.bioconductor.org/packages/release/bioc/manuals/dada2/man/dada2.pdf) for descriptions of each function. 

>For example, the manual explains:
>The distribution of quality scores at each position is shown as a grey-scale heat map, with dark
>colors corresponding to higher frequency. The plotted lines show positional summary statistics:
>green is the mean, orange is the median, and the dashed orange lines are the 25th and 75th quantiles.
>If the sequences vary in length, a red line will be plotted showing the percentage of reads that extend
>to at least that position

Now we will filter and trim. First we will make lists of the outputs files we will produce: `filtFs` and `filtRs`.

```R
# list filtered files in subdirectory: called filtered

# make a list of filtered output files for forward reads:
filtFs <- file.path(outpath, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) # forward reads
names(filtFs) <- sample.names

# make a list of filtered output files for reverse reads:
filtRs <- file.path(outpath, "filtered", paste0(sample.names, "_R_filt.fastq.gz")) # reverse reads
names(filtRs) <- sample.names

````

The `filterAndTrim` function

> Parameters to think about are:
>`trimleft`
>`trunclen`
>`maxN`
>`maxEE`
>`truncQ`


``` R
# Action required: Change trunclen=c(xx,xx) to match what you want to truncate your forward and reverse reads to. Similarly edit trimLeft=c(xx,xx), maxEE=x and truncQ=x. ##change_me
# https://rdrr.io/bioc/dada2/man/filterAndTrim.html
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(10,10), truncLen=c(200,160),
              maxN=0, maxEE=c(2,2), truncQ=c(2,2), rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=TRUE
head(out) # check how many reads have been lost in filtering step
# save out file in case environment shuts down
write.table(out,paste(outpath,"/out_save.tsv",sep=""),sep="\t",row.names=T)
```

At this point we can look at the summary of the trimming step by looking at the varable `out`. To do this type `out` into the terminal.
>Or, for ease, print "out" ordered by column two, which is the number of reads left: `out[order(out[,2], decreasing=T),]`

Questions: Have many reads been lost? Are there any samples with too few reads to proceed with?

Next we can make more plots of the read files, this time after they have been trimmed.

``` R
# make a pdf of the read quality output plot_read_quality_F.pdf
pdf(file = paste(outpath,"/plot_read_quality_filtered_F.pdf",sep=""), 7, 7)
plotQualityProfile(filtFs[1:12])
dev.off()

# make a pdf of the read quality output plot_read_quality_R.pdf
pdf(file = paste(outpath,"/plot_read_quality_filtered_R.pdf",sep=""), 7, 7)
plotQualityProfile(filtRs[1:12])
dev.off()
```

Sometimes files lose all reads at this step. If that happens you need to exclude them from the rest the analysis.

Or, if there are samples with some reads, but you would like to exclude them from the rest of the analysis, you can remove them with the code below.

More specifically, the variable `removed_samples` is made first, and is populated by a list of samples with zero reads left (`which(as.data.frame(out)[,2] ==0`). To add samples to this list you need to use `c(removed samples, "sample_file_1.fastq.gz")`
e.g.:
`removed_samples<- c(removed_samples, "P1a_1.fastq.gz", "P2a_1.fastq.gz", "P1e_1.fastq.gz", "P1c_1.fastq.gz")`

``` R
# if statement checks if samples need to be removed and following statements remove them
removed_samples<-row.names(as.data.frame(out)[(which(as.data.frame(out)[,2] ==0)),])

# modify this line to remove samples manually.
# here, sample Ts29 has fewer reads than the rest of the dataset so we will remove it
# note that in the table "out" full file names with extensions are required e.g:"Ts29_1.fastq.gz"
removed_samples<-c(removed_samples,"Ts29_1.fastq.gz")

# clean up file names
# change this if your file name is not structured as SAMPLENAME_fileextension. If there are "_" in your samplename you will need to change this.
removed_samples.names <- sapply(strsplit(basename(removed_samples), "_"), `[`, 1) # check this works for your data

# remove samples from filtFs, filtRs, sample.names
if(length(removed_samples) != 0){
  filtFs<-filtFs[(!(names(filtFs) %in% removed_samples.names))] # remove from list of forward filtered read files
  filtRs<-filtRs[(!(names(filtRs) %in% removed_samples.names))] # remove from list of reverse filtered read files
  sample.names<-sample.names[!(sample.names) %in% removed_samples.names] # removes from sample names list
  out<-out[which(!(row.names(out)%in% removed_samples)),] # removes from out
}

```

Stop here and check the filtering and trimming step before moving forward!
Do you need to repeat with new parameters?

So, if it seems we have lost too many reads, or not enough, or the quality scores are still too low, we need to read try again with new parameters.

We can check the read quality using FastQC again to be sure:

### 4.3.2 Run [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [Multiqc](https://www.bing.com/search?q=multiqc&cvid=41c889025cf2474eae8615972425b160&aqs=edge..69i57j0l8j69i11004.1847j0j4&FORM=ANAB01&PC=U531) [optional]

#### Open a new command-line session:

So that we don't interrupt the R session we have been working in, let's start a new session.

Option 1 - You can exit your screen using `ctrl + a then d`, and start a new one `screen -S fastqc`.
Then to return you can use `ctrl + a then d`, and re-enter the `meta` screen using `screen -r meta`

Or Option 2 - You can start a new session using Mobaxterm and work there

```bash
# start another interactive session in a different session
	msub -I -q S30 -l procs=8,walltime=0:30:00

# Activate conda environment:
	conda activate fastqc

# Set paths
	indir=~/metabarcoding_ws/outputs/dada2/filtered
	outdir=~/metabarcoding_ws/outputs/dada2/

# Make a directory to store fastqc files
	mkdir $outdir/filtered_fastqc_out

# Run fastqc:
	fastqc $indir/*.fastq.gz -t 8 -o $outdir/filtered_fastqc_out
# runs fastqc on all fastq.gz n current directory. -t threads to use

# Run multiqc on fastqc output folder.
#(It automatically detects fastqc outputs)
# change `--title` if you wish
	multiqc $outdir/filtered_fastqc_out/* -o $outdir/filtered_fastqc_out/multiqc_fastqc --interactive --title fastqc_filtered_interactive
```

If you are happy to proceed, go back to the R session and start to calculate error rates:

## 4.4 Error rates
### Back in R:

Dada2 relies on error models calculated for every dataset. the function `learnErrors()` will calculate these and store them in `errF` and `errR`

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

The function `derepFastq()` will group all the identical reads

``` R
derep_forward <- derepFastq(filtFs, verbose=TRUE)
# name the derep-class objects by the sample names
names(derep_forward) <- sample.names

derep_reverse <- derepFastq(filtRs, verbose=TRUE)
# name the derep-class objects by the sample names
names(derep_reverse) <- sample.names
```

## 4.6 Sample Inference/ Denoise reads

This utilises the main function of Dada2 `dada()`

Here, dada2 will use error models calculated for forward and reverse reads separaetly to determine which ASV's are likely to be real biological sequences or if they have occured due to errors. This is done separately since assessing overlapped regions between paired reads is less accurate than analysing each read in the pair directly. 

After denoising, ASV's are then merged by aligning paired reads. 

Then, an ASV table is created `seqtab` using function `makeSequenceTable`. 

We can check that the merged reads are the length we expected using `getSequences()`.

``` R
# Forward - this clusters forward reads into ASV's. later you will merge forward and reverse reads
dadaFs <- dada(derep_forward, err=errF, multithread=TRUE)
# Reverse - this clusters reverse reads into ASV's.
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
## 4.7 Remove chimeras

We need to remove chimeras. A function in dada2 called `removeBimeraDenovo` will do this by removing reads that look like a hybrid of two others in the dataset. Chimeras made from more than two sequences merging are rare.

``` R
# run removeBimeraDenovo
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# check size of data.frame
dim(seqtab.nochim)

# this line will print a message using the outputs:
print(paste(round(sum(seqtab.nochim)/sum(seqtab),4)*100,"% of ASV's remain after chimera removal. Remaining ASV's:",dim(seqtab.nochim[2]),sep=""))
```

## 4.8 Removal of negative controls (Optional)

If your negative control doesn't contain reads, or you have not included a negative control in your study, skip this step. 

If you have included negative controls in your experiment, and they have enough reads remaining at this stage that ASV's have been found in those samples, you can remove those ASV's from the rest of the dataset. This is assuming that these ASV's have been introduced to the rest of your data through contamination. This method may be too simple if there are contaminants that are also present in the sample.

You will need to make a list of negative control samples e.g: `negative_controls<-c("T16s", "T20s","T21s")`. You will only need the sample names for this stage, not the file name as before when we removed samples.

Here we will pretend that sample `Ts16s` is a negative control.

``` R
# make a list of the negative control samples
negative_controls<-c("T16s")

# save a new version of the ASV table `seqtab.nochim`
seqtab.nochim.c<-seqtab.nochim

# this loop will loop through each negative control (f) in the list (negative controls) and remove all ASV's found in that sample (to_remove) from the rest of the samples
for(f in negative_controls){
	# print which negative control is being used in the loop righ now
  print(f)
	# make a list of rows to drop so we can count how many
  to_remove<-dim(seqtab.nochim.c[,which(seqtab.nochim.c[which(row.names(seqtab.nochim.c)==f),]>0)])[2]
	# keeps rows where negative control is zero
  seqtab.nochim.c<-seqtab.nochim.c[,which(seqtab.nochim.c[which(row.names(seqtab.nochim.c)==f),]<=0)]
	# print a message to say how many ASV's were dropped and kept
  print(paste(dim(seqtab.nochim.c)[2]," ASV's remain after removing", to_remove," ASV's present in negative controls:", f))
}

# if we would like to keep this new ASV table, we can replace the old dataframe seqtab.nochim with the new one seqtab.nochim.c and proceed:
seqtab.nochim<-seqtab.nochim.c
rm(seqtab.nochim.c)

# we won't do this here because we didn't use a true negative control

```

Now that we have removed chimeras and optionally removed negative controls we can save our final ASV table which is saved as `seqtab.nochim`.

``` R
#sample.names <- 

sapply(strsplit(seqtab.nochim, "_"), `[`, 1)

names(seqtab.nochim)

# save result
write.table(seqtab.nochim, file = paste(outpath,"/asv_table.dada2.tsv",sep=""),sep="\t",row.names=T)

# Save seqtab.nochim as an RDS file - this will make reading the table back into R much easier
saveRDS(seqtab.nochim, paste(outpath,"/seqtab_nochim.rds",sep=""))

# save csv for MicrobiomeAnalyst
# transpose
seqtab.nochim.df <- t(seqtab.nochim)
# make a column for the samples names called "#NAMES"
seqtab.nochim.df <- cbind(row.names(seqtab.nochim.df),seqtab.nochim.df)

# set first column name as "#NAMES" for MicrobiomeAnalyst
colnames(seqtab.nochim.df)[1]<-"#NAMES"

# write output
write.csv(seqtab.nochim.df, file = paste(outpath,"/asv_table.dada2.ma.csv",sep=""),row.names=F)
```

## 4.9 Track data across each step

Let's build a table that shows how many reads/ASV's are left after each step.

``` R
# set getN function - this function will count for us
getN <- function(x) sum(getUniques(x))

# Combine filtering, dada2, and chimera removal data into one dataframe
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), rowSums(seqtab.nochim.c))

# set some column names
colnames(track) <- c("no_reads_input", "no_reads_after_filter", "no_ASV_denoisedF", "no_ASV_denoisedR", "no_ASV_merged", "no_ASV_nonchim","no_ASV_neg_cntrl_removed")
rownames(track) <- sample.names

# check output
head(track) # head shows the first 10 entries

# write to file:
write.table(track, paste(outpath,"/track_reads.tsv", sep=""), sep="\t")
```

>We can view this table in order of reads remain after negative controls are removed as follows:
`track[order(track[,7], decreasing=T),]`

> Question: How many ASV's survived the denoising and merging processes?

## 4.10 Taxonomy

Next, we will assign taxa to each ASV we have detected. We will compare our ASV's to he reference database SILVA.

Functions for assigning taxonomy are `taxa_with_bootstraps()` and `taxa()`. Here we use `taxa_with_bootstraps()`.

``` R
# dataset - SILVA:
taxa_with_bootstraps <- assignTaxonomy(seqtab.nochim, "~/metabarcoding_ws/db/silva_nr99_v138.1_train_set.fa.gz", outputBootstraps=TRUE, multithread=TRUE)

# we can save the taxa dataframe without the bootstraps:
taxa_all <- taxa_with_bootstraps$tax

# optionally add species:
taxa_species <- addSpecies(taxa_with_bootstraps$tax, "~/metabarcoding_ws/db/silva_species_assignment_v138.1.fa.gz")

# write to file
write.csv(taxa_with_bootstraps, file = paste(outpath,"/taxa_with_bootstraps.csv",sep=""),row.names=T)
write.csv(taxa_all, file = paste(outpath,"/taxa.csv",sep=""),row.names=T)
write.csv(taxa_species, file = paste(outpath,"/taxa_species.csv",sep=""),row.names=T)

# write for MicrobiomeAnalyst
taxa_df <- cbind(row.names(taxa_all),taxa_all)
colnames(taxa_df)[1]<-"#TAXONOMY"

# write - important that "row.names=F" this time
write.csv(taxa_df, file = paste(outpath,"/taxa.ma.csv",sep=""),row.names=F)
```

**Now you have completed the Dada2 aspect of the analysis!**

>You should now  have been left with output files in your `outpath`: `/home/$USER/metabarcoding_ws/outputs/dada2/`
>
>**Most important:**
>
>ASV file: `asv_table.dada2.tsv`
>
>Taxonomy files: `taxa.csv`, `taxa_with_bootstraps.csv`, and `taxa_species.csv`
>
>MicrobiomeAnalyst format files :  `asv_table_dada2.ma.csv`, `taxa.ma.csv`
>
>**Additional:**
>
>Filtered reads directory: `filtered/`
>
>Read quality pdfs (pre filtering): `plot_read_quality_F.pdf`, `plot_read_quality_R.pdf`
>
>Read quality pdfs (post filtering): `plot_read_quality_filtered_F.pdf`, `plot_read_quality_filtered_R.pdf`
>
>R formatted ASV table: `seqtab_nochim.rds`
>
>Track reads file: `track_reads.tsv`
>
>Saved out file (can be discarded now): `out_save.tsv`
>
>Error rate plot pdfs: `plot_errors_F.pdf`, `plot_errors_R.pdf`

## 5. Alternatively run [dadaist2](https://github.com/quadram-institute-bioscience/dadaist2) (Optional)

This is a command line wrapper for dada2.

#### In the command-line:
``` bash
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
```

# 6. Visualisation and statistics:

## 6.1 [MicrobiomeAnalyst](https://www.microbiomeanalyst.ca/)

MicrobiomeAnalyst is a useful tool for exploring your data (although it can be tricky to use if your data isn't in the correct format).

First download the metadata file:

```
wget https://raw.githubusercontent.com/nmc97/metabarcoding_at_cefas_tutorial/main/scripts/workshop_1/data/metadata.csv
```

Then use the files we produced for MicrobiomeAnalyst:  `asv_table_dada2.ma.csv`, `taxa.ma.csv`

Input these data here:
https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/upload/OtuUploadView.xhtml

**Dadaist2**

To submit the files in the Dadaist2 MicrobiomeAnalyst directory you may first need to edit some of them.

The table.txt file needs to have file extensions removed from the sample names (use find and replace).

The taxa file doesn't contain enough delimiters. Open in excel and save as csv file.

## 6.2 Phyloseq

For this you need R. It would be best to set up R studio on your computer and run this there.

You will need to create a Phyloseq object from the outputs of Dada2 or Dadaist2

## 6.3 Microbiome R package

For this you need R. It would be best to set up r studio on your computer and run this there.

You will need to create a Phyloseq object from the outputs of Dada2 or Dadaist2

https://microbiome.github.io/tutorials/
