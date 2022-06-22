##=================================================================================#
##
## Rscript to run dada2
##
## Expects demultiplexed reads stored in a single folder (path)
##
## Author: Nicola Coyle, Diana Minardi, David Ryder?
## Institution: Cefas
## Contact: nicola.coyle@cefas.co.uk
##
## Adapted from dada2 tutorial: https://benjjneb.github.io/dada2/tutorial.html
##
##=================================================================================#

#=========================#
# Prerequisites
# - DADA2
#   install DADA2 - https://metabarcoding-at-cefas-tutorial.readthedocs.io/en/latest/dada2.html#dada2-background-installation
# - dplyr
#=========================#

#=========================#
# load libraries
#=========================#

library(dada2)
library(dplyr)

#=========================#
# Set paths and parameters here
#=========================#
# set outpath
outpath <- "/home/nc07/projects/metabarcoding/programs/dada2/dadaist2/example_out_dada2script"
# set path to data directory
path <- "/home/nc07/projects/metabarcoding/programs/dada2/dadaist2/data/16S_rename" # trimmed reads

# set file extension
ext_R1 <- "_R1.fastq.gz"
ext_R2 <- "_R2.fastq.gz"

# set taxonomy reference dataset file
tax_file <- "/home/nc07/projects/metabarcoding/programs/dada2/silva_nr99_v138.1_train_set.fa.gz"
tax_file_species <- "/home/nc07/projects/metabarcoding/programs/dada2/silva_species_assignment_v138.1.fa.gz" # optional
dataset_type <- "silva" # type silva or pr2

# use forward reads only - set to TRUE (sometimes the reverse reads are too poor to use)
forward <- FALSE

# run filter step TRUE/ FALSE (you may have already filtered and trimmed your reads)
filter_reads <- TRUE

# set filter parameters
trimLeft <- c(17,17) # change to a single number if only using forward reads
truncLen <- c(280,280) # change to a single number if only using forward reads
maxN <- 0
maxEE <- 2
truncQ <- 2
rm.phix <- TRUE
compress <- TRUE
multithread <- TRUE  # On Windows set multithread=FALSE

MicrobiomeAnalyst <- TRUE

#=========================#
# setup directories and files
#=========================#

# if statement makes output folder if doesn't exist
if (file.exists(outpath)) {
 cat("Output directory already exists")
} else {
 dir.create(outpath)
 cat("Output directory created")
}

# check files exist in read file directory
if(length(list.files(path))==0){
  print("No files in directory provided.")
} else {
  print(paste("directory contains",length(list.files(path)),"files",sep=" "),)
  print("First 10:")
  print(head(list.files(path)))
}

# Using the file extension pattern given above "ext_1" and "ext_R2", read files are found
fnFs <- sort(list.files(path, pattern=ext_R1, full.names = TRUE))
if(forward == FALSE){fnRs <- sort(list.files(path, pattern=ext_R2, full.names = TRUE))}
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq.gz
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) # check this works for your data!!
if(length(fnFs)==length(fnRs)){
  print(paste(length(sample.names)," samples found",sep=""))
} else {
  if(forward == FALSE){
    print("Unequal number of paired reads.")
  } else {
    print(paste(length(sample.names)," samples found",sep=""))
  }
}

#=========================#
# Filter and trim
#=========================#

# if you set the filter reads parameter to be TRUE reads will be Filtered
# please check the filter parameters

if(filter_reads ==TRUE){
  # Place filtered files in filtered/ subdirectory
  filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) # forward reads
  names(filtFs) <- sample.names
  if(forward == FALSE){filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz")) # reverse reads
  names(filtRs) <- sample.names}

  # check filter and trim parameters !
  # this function will create a directory called "filtered" in the read directory
  if(forward == FALSE){
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=trimLeft, truncLen=truncLen, maxN=maxN, maxEE=maxEE, truncQ=truncQ, rm.phix=rm.phix, compress=compress, multithread=multithread)
  } else {
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=trimLeft, truncLen=truncLen, maxN=maxN, maxEE=maxEE, truncQ=truncQ, rm.phix=rm.phix, compress=compress, multithread=multithread)
  }

  print(paste("Filtered reads can be found in directory: ",path, "/filtered",sep=""))

  print("Review percentage of lost reads:")
  print(as.data.frame(out) %>% mutate(pc.reads.lost=(((reads.in-reads.out)/reads.in)*100)))

  # save out file in case environment shuts down
  write.table(out,paste(outpath,"/out_save.tsv",sep=""),sep="\t",row.names=T)

  # Sometimes files loose all reads at this step. If that happens you need to exclude them from the rest the analysis.
  # if statement checks if samples need to be removed and following statements remove them
  removed_samples<-row.names(as.data.frame(out)[(which(as.data.frame(out)[,2] ==0)),])
  removed_samples.names <- sapply(strsplit(basename(removed_samples), "_"), `[`, 1) # check this works for your data
  if(length(removed_samples) != 0){
    filtFs<-filtFs[(!(names(filtFs) %in% removed_samples.names))] # remove from list of forward filtered read files
    if(forward == FALSE){filtRs<-filtRs[(!(names(filtRs) %in% removed_samples.names))]} # remove from list of reverse filtered read files
    sample.names<-sample.names[!(sample.names) %in% removed_samples.names]
    out<-out[which(!(row.names(out)%in% removed_samples)),]
  }
} else {
  print("Skipping filter step.")
  filtFs<-fnFs
  if(forward == FALSE){filtRs<-fnRs}
}

#=========================#
# Error rates
#=========================#

# calculate error rates for forward and reverse reads
errF <- learnErrors(filtFs, multithread=TRUE)
if(forward == FALSE){errR <- learnErrors(filtRs, multithread=TRUE)}

# make a pdf of the error output plot_errors.pdf
pdf(file = paste(outpath,"/plot_errors_F.pdf",sep=""), 7, 7)
plotErrors(errF, nominalQ=TRUE)
dev.off()
if(forward == FALSE){
pdf(file = paste(outpath,"/plot_errors_R.pdf",sep=""), 7, 7)
plotErrors(errR, nominalQ=TRUE)
dev.off()
}

#=========================#
# Dereplicate reads
#=========================#

derep_forward <- derepFastq(filtFs, verbose=TRUE)
# name the derep-class objects by the sample names
names(derep_forward) <- sample.names

if(forward == FALSE){
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
# name the derep-class objects by the sample names
names(derep_reverse) <- sample.names
}

#=========================#
# Sample inference
#=========================#

# Forward
dadaFs <- dada(derep_forward, err=errF, multithread=TRUE)
# Reverse
if(forward == FALSE){dadaRs <- dada(derep_reverse, err=errR, multithread=TRUE)}

# merged paired read data
if(forward == FALSE){
  clusters <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
} else {
  clusters<-dadaFs
}
# Inspect the merger data.frame from the first sample
head(clusters[[1]])

# amplicon sequence variant table (ASV)
seqtab <- makeSequenceTable(clusters)
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
# track reads
#=========================#

# use this to determine which samples would be worth keeping in merged run

# set getN function
getN <- function(x) sum(getUniques(x))

# Combine filtering, dada2, and chimera removal data into one dataframe
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
if(forward == FALSE){
  track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(clusters, getN), rowSums(seqtab.nochim))
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sample.names
} else {
  track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
  colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
  rownames(track) <- sample.names
}

print("Top of read tracking table:")
print(head(track)) # check output

# write to file:
write.table(track, paste(outpath,"/track_reads.tsv", sep=""), sep="\t")
print(paste("written to file ", paste(outpath,"/track_reads.tsv", sep=""),sep=""))

#=========================#
# Taxonomy
#=========================#

# dataset - if PR2 database used adding taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"):
if(dataset_type == "silva"){
  print(paste("Assigning taxa using a SILVA database: ",tax_file, sep=''))
  taxa_with_bootstraps <- assignTaxonomy(seqtab.nochim, tax_file, outputBootstraps=TRUE, multithread=TRUE)
  taxa_all <- taxa_with_bootstraps$tax
    # optionally add species:
    print(paste("Assigning species using a SILVA database: ",tax_file_species, sep=''))
    if(tax_file_species != ''){
    taxa_species <- addSpecies(taxa_with_bootstraps[[1]], tax_file_species)
    }
} else if (dataset_type == "pr2"){
  print(paste("Assigning taxa using a PR2 database: ",tax_file, sep=''))
  taxa_with_bootstraps <- assignTaxonomy(seqtab.nochim, tax_file, taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"), outputBootstraps=TRUE, multithread=TRUE) ##change_me
  taxa_all <- taxa_with_bootstraps$tax
} else {
  print(paste("Assigning taxa using an unknown database type: ",tax_file, sep=''))
  taxa_with_bootstraps <- assignTaxonomy(seqtab.nochim, tax_file, outputBootstraps=TRUE, multithread=TRUE)
  taxa_all <- taxa_with_bootstraps$tax
}

# write to files
write.table(taxa_with_bootstraps, file = paste(outpath,"/taxa_with_bootstraps.tsv",sep=""),sep="\t",row.names=T)
print(paste("Taxa bootstrap table written to file ", paste(outpath,"/taxa_with_bootstraps.tsv", sep=""),sep=""))
write.table(taxa_all, file = paste(outpath,"/taxa.tsv",sep=""),sep="\t",row.names=T)
print(paste("Taxa table written to file ", paste(outpath,"/taxa.tsv", sep=""),sep=""))
# Only write to a species file if a species reference dataset was provided:
if(tax_file_species != ''){write.table(taxa_species, file = paste(outpath,"/taxa.species.tsv",sep=""),sep="\t",row.names=T)}
print(paste("Taxa species table written to file ", paste(outpath,"/taxa.species.tsv", sep=""),sep=""))

# create files to go into MicrobiomeAnalyst
if(MicrobiomeAnalyst==TRUE){
  MA_outpath<-paste(outpath,"/MicrobiomeAnalyst",sep="")
  if (file.exists(MA_outpath)) {
  cat("Output directory already exists")
  } else {
  dir.create(MA_outpath)
  cat("Output directory created")
  }

  ma_seqtab.nochim<- tibble::rownames_to_column(as.data.frame(t(seqtab.nochim)), "#NAME")
  write.table(ma_seqtab.nochim, file = paste(MA_outpath,"/asv_table.MA.dada2.csv",sep=""),sep=",",row.names=F, quote=F)

  ma_taxa_all<- tibble::rownames_to_column(as.data.frame((taxa_all)), "#TAXONOMY")
  write.table(ma_taxa_all, file = paste(MA_outpath,"/taxa.MA.csv",sep=""),sep=",",row.names=F, quote=F, na = "")
  print(paste("Taxa table written to file ", paste(MA_outpath,"/taxa.MA.csv", sep=""),sep=""))
  # Only write to a species file if a species reference dataset was provided:
  if(tax_file_species != ''){
    ma_taxa_species<- tibble::rownames_to_column(as.data.frame((taxa_species)), "#TAXONOMY")
    write.table(ma_taxa_species, file = paste(MA_outpath,"/taxa.species.MA.csv",sep=""),sep=",",row.names=F, quote=F, na = "")
    print(paste("Taxa species table written to file ", paste(MA_outpath,"/taxa.species.MA.tsv", sep=""),sep=""))
  }
}

print(done)
