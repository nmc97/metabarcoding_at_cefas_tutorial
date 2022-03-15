##=================================================================================#
##
## Rscript to run dada2
##
## Expects demultiplexed reads stored in a single folder
##
## Author: Nicola Coyle, Diana Minardi, David ryder?
## Institution: Cefas
## Contact: nicola.coyle@cefas.co.uk
##
## Adapted from dada2 tutorial: https://benjjneb.github.io/dada2/tutorial.html
##
##=================================================================================#

# note: search for string ##change_me to find lines containing things you need to change for your data

#=========================#
# load libraries
#=========================#

library(dada2)
library(dplyr)
library(reshape2)
library(tibble) # optional

#=========================#
# setup file paths
#=========================#

# set outpath
outpath <- "/path/to/output/directory" ##change_me
# set path to data directory and list files
path <- "/path/to/read/directory" # trimmed reads ##change_me

# if statement makes output folder if doesn't exist
if (file.exists(outpath)) {
 cat("The folder already exists")
} else {
 dir.create(outpath)
}

# check files exist in read file directory
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE)) ##change_me if the file names are different to this structure
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE)) ##change_me if the file names are different to this structure
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#=========================#
# Filter and trim
#=========================#

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) # forward reads
names(filtFs) <- sample.names
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz")) # reverse reads
names(filtRs) <- sample.names

# Action required: Change trunclen=c(xx,xx) to match what you want to truncate your forward and reverse reads to. ##change_me
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(17,17), truncLen=c(280,280),
              maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out) # check how many reads have been lost in filtering step
# save out file in case environment shuts down
write.table(out,paste(outpath,"/out_save_F1.tsv",sep=""),sep="\t",row.names=F)


# Sometimes files loose all reads at this step. If that happens you need to exclude them from the rest the analysis. Here is an example of how you could do so with two samples "sample_name1" and "sample_name2"
# ##change_me : uncomment the following lines and change samples names to match your data
# filtFs<-filtFs %>% filter(filtFs!="sample_name1" & filtFs!="sample_name2") # remove from list of forward filtered read files
# filtRs<-filtRs %>% filter(filtRs!="sample_name1" & filtRs!="sample_name2") # remove from list of reverse filtered read files
# sample.names<-sample.names[(sample.names)!="sample_name1" $ (sample.names)!="sample_name2"] # remove from sample name list

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

# Forward
dadaFs <- dada(derep_forward, err=errF, multithread=TRUE)
# Reverse
dadaFs <- dada(derep_reverse, err=errF, multithread=TRUE)

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

write.table(seqtab.nochim, file = paste(outpath,"/taxa_with_bootstraps.tsv",sep=""),sep="\t",row.names=F)

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
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
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
# taxa_with_bootstraps <- assignTaxonomy(seqtab.nochim, "/path/to/ref/dataset/silva_nr99_v138.1_train_set.fa.gz", outputBootstraps=TRUE, multithread=TRUE) ##change_me
# taxa_all <- taxa_with_bootstraps$taxa

# dataset - Pr2:
# must define taxlevels for Pr2
taxa_with_bootstraps <- assignTaxonomy(seqtab.nochim, "/path/to/ref/dataset/pr2_version_4.14.0_SSU_dada2.fasta.gz", taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"), outputBootstraps=TRUE, multithread=TRUE) ##change_me
taxa_all <- taxa_with_bootstraps$taxa

# write to files
write.table(taxa_with_bootstraps, file = paste(outpath,"/taxa_with_bootstraps.tsv",sep=""),sep="\t",row.names=F)
write.table(taxa, file = paste(outpath,"/taxa.tsv",sep=""),sep="\t",row.names=F)
