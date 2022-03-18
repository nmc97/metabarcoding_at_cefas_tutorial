#! bash
##=================================================================================#
##
## Rscript to run quality testing of read files
##
## Expects a folder containing gzipped read files
##
## Author: Nicola Coyle
## Date: March 2022
## Institution: Cefas
## Contact: nicola.coyle@cefas.co.uk
##
##=================================================================================#

# uses fastqc and multiqc so create a mamba environment with these programs
# mamba create -n fastqc fastqc MultiQC
# check they both work and them proceed

# start interactive session (alterantively run a job)
msub -I -q S30 -l procs=12,walltime=3:00:00 # change depending on what yo think you need.

conda activate fastqc

cd /path/to/directory/containing/read/files ##change_me

fastqc *.fastq.gz -t 12 # runs fastqc on all fastq.gz n current directory
mkdir fastqc_out
mv *fastqc.* fastqc_out/ # move output files

# run multiqc on fastqc output folder. it automatically detects fastqc outputs
multiqc fastqc_out/* -o ./fastqc_out/multiq_fastqc --title fastqc_flat # change title
# force interative if there are lots of files
multiqc fastqc_out/* -o ./fastqc_out/multiqc_fastqc --interactive --title fastqc_interactive

# download and inspect output files
