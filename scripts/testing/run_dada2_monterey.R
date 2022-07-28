##=================================================================================#
##
## Rscript to run dada2
##
## Expects demultiplexed reads stored in a single folder
##
## Author: Nicola Coyle, Diana Minardi, David Ryder
## Institution: Cefas
## Contact: nicola.coyle@cefas.co.uk
##
## Adapted from dada2 tutorial: https://benjjneb.github.io/dada2/tutorial.html
##
## using data from PRJNA726280
##=================================================================================#

# before using:
# install dada2 - https://metabarcoding-at-cefas-tutorial.readthedocs.io/en/latest/dada2.html#dada2-background-installation

# note: search for string ##change_me to find lines containing things you need to change for your data

# /home/nc07/projects/metabarcoding/workshop/monterey

# In bash:
```
conda activate fastqc

indir=/home/nc07/projects/metabarcoding/workshop/monterey/16S/fastq
outdir=/home/nc07/projects/metabarcoding/workshop/monterey/outputs/16S/

# move to directory containing read files ##change_me
cd /path/to/directory/containing/read/files ##change_me
mkdir /home/nc07/projects/metabarcoding/workshop/monterey/outputs/16S/fastqc_out # make directory to store outputs

fastqc /home/nc07/projects/metabarcoding/workshop/monterey/16S/fastq/*.fastq.gz -t 12 -o /home/nc07/projects/metabarcoding/workshop/monterey/outputs/16S/fastqc_out# runs fastqc on all fastq.gz n current directory. -t threads to use

# run multiqc on fastqc output folder. it automatically detects fastqc outputs
multiqc /home/nc07/projects/metabarcoding/workshop/monterey/outputs/16S/fastqc_out/* -o /home/nc07/projects/metabarcoding/workshop/monterey/outputs/16S/fastqc_out/multiq_fastqc --title fastqc # change title
# force interative if there are lots of files (it will tell you if it wrote flat files instead)
# multiqc /home/nc07/projects/metabarcoding/workshop/monterey/outputs/16S//* -o /home/nc07/projects/metabarcoding/workshop/monterey/outputs/16S/fastqc_out/multiqc_fastqc --interactive --title fastqc_interactive
```

#=========================#
# load libraries
#=========================#

library(dada2)

#=========================#
# setup file paths
#=========================#

# set outpath
outpath <- "/home/nc07/projects/metabarcoding/workshop/monterey/outputs/16S" ##change_me - this is where your output files will go and it will be created for you if the next directory above exists
# set path to data directory and list files
path <- "/home/nc07/projects/metabarcoding/workshop/monterey/16S/fastq" # trimmed reads ##change_me

# if statement makes output folder if doesn't exist
if (file.exists(outpath)) {
 cat("The folder already exists")
} else {
 dir.create(outpath)
}

# check files exist in read file directory
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE)) ##change_me if the file names are different to this structure
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE)) ##change_me if the file names are different to this structure
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#=========================#
# Filter and trim
#=========================#

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(outpath, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) # forward reads
names(filtFs) <- sample.names
filtRs <- file.path(outpath, "filtered", paste0(sample.names, "_R_filt.fastq.gz")) # reverse reads
names(filtRs) <- sample.names

# Action required: Change trunclen=c(xx,xx) to match what you want to truncate your forward and reverse reads to. Similarly edit trimLeft=c(xx,xx), maxEE=x and truncQ=x. ##change_me
# https://rdrr.io/bioc/dada2/man/filterAndTrim.html
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(10,10), truncLen=c(150,150),
              maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=TRUE
head(out) # check how many reads have been lost in filtering step
# save out file in case environment shuts down
write.table(out,paste(outpath,"/out_save_F1.tsv",sep=""),sep="\t",row.names=T)

# Sometimes files loose all reads at this step. If that happens you need to exclude them from the rest the analysis.
# if statement checks if samples need to be removed and following statements remove them
removed_samples<-row.names(as.data.frame(out)[(which(as.data.frame(out)[,2] ==0)),])
removed_samples<-c(removed_samples,"ESPCTD16S47_1.fastq.gz","ESPCTD16S80_1.fastq.gz","ESPCTD16S42_1.fastq.gz","ESPCTD16S36_1.fastq.gz")
removed_samples.names <- sapply(strsplit(basename(removed_samples), "_"), `[`, 1) # check this works for your data
if(length(removed_samples) != 0){
  filtFs<-filtFs[(!(names(filtFs) %in% removed_samples.names))] # remove from list of forward filtered read files
  filtRs<-filtRs[(!(names(filtRs) %in% removed_samples.names))] # remove from list of reverse filtered read files
  sample.names<-sample.names[!(sample.names) %in% removed_samples.names]
  out<-out[which(!(row.names(out)%in% removed_samples)),]
}

# stop here and check the filtering step before moving forward! Repeat with new parameters if needed.

#=========================#
# Error rates
#=========================#

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

#=========================#
# Dereplicate reads
#=========================#

derep_forward <- derepFastq(filtFs, verbose=TRUE)
# name the derep-class objects by the sample names
names(derep_forward) <- sample.names

derep_reverse <- derepFastq(filtRs, verbose=TRUE)
# name the derep-class objects by the sample names
names(derep_reverse) <- sample.names

#=========================#
# Sample inference
#=========================#

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

#=========================#
# remove chimeras
#=========================#

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) # check size of data.frame
sum(seqtab.nochim)/sum(seqtab) # check proportion of reads are left after chimeras removed

write.table(seqtab.nochim, file = paste(outpath,"/asv_table.dada2.tsv",sep=""),sep="\t",row.names=T)

#=========================#
# remove negative control ASVs
#=========================#

negative_controls<-c("SRR14366259","ESPCTD16S93")
seqtab.nochim.c<-seqtab.nochim
for(f in negative_controls){
  print(f)
  to_remove<-dim(seqtab.nochim.c[,which(seqtab.nochim.c[row.names(seqtab.nochim.c)==f,]>0)])[2]
  seqtab.nochim.c<-seqtab.nochim.c[,which(seqtab.nochim.c[row.names(seqtab.nochim.c)==f,]<=0)]
  print(paste(dim(seqtab.nochim.c)[2]," ASV's remain after removing", to_remove," ASV's present in negative controls:", f))
}
# seqtab.nochim<-seqtab.nochim.c
# rm(seqtab.nochim.c)
# read.table(paste(outpath,"/asv_table.dada2.tsv",sep=""),sep="\t",header=T)

#=========================#
# track reads
#=========================#

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

#=========================#
# Taxonomy
#=========================#

# decide what reference dataset to use
# examples:

# dataset - SILVA:
taxa_with_bootstraps <- assignTaxonomy(seqtab.nochim, "~/projects/metabarcoding/programs/dada2/silva_nr99_v138.1_train_set.fa.gz", outputBootstraps=TRUE, multithread=TRUE) ##change_me
taxa_all <- taxa_with_bootstraps$tax
# optionally add species:
taxa_species <- addSpecies(taxa_with_bootstraps$tax, "~/projects/metabarcoding/programs/dada2/silva_nr99_v138.1_train_set.fa.gz/silva_species_assignment_v132.fa.gz")

# dataset - Pr2:
# must define taxlevels for Pr2
#taxa_with_bootstraps <- assignTaxonomy(seqtab.nochim.c, "/path/to/ref/dataset/pr2_version_4.14.0_SSU_dada2.fasta.gz", taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"), outputBootstraps=TRUE, multithread=TRUE) ##change_me
#taxa_all <- taxa_with_bootstraps$tax

# write to files
write.table(taxa_with_bootstraps, file = paste(outpath,"/taxa_with_bootstraps.tsv",sep=""),sep="\t",row.names=T)
write.table(taxa_all, file = paste(outpath,"/taxa.tsv",sep=""),sep="\t",row.names=T)
write.table(taxa_species, file = paste(outpath,"/taxa_species.tsv",sep=""),sep="\t",row.names=T)
