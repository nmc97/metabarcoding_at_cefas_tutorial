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
# create a conda environment with: `mamba create -n fastqc fastqc MultiQC`
# check both programs work and then proceed

# start interactive session (alterantively run a job)
msub -I -q S30 -l procs=12,walltime=2:00:00 # change depending on what yo think you need.

conda activate fastqc

#'cd /path/to/directory/containing/read/files ##change_me
cd /home/nc07/projects/metabarcoding/programs/dada2/dadaist2/data/16S_rename/ ##change_me
mkdir fastqc_out # make directory to store oututs

fastqc *.fastq.gz -t 12 -o fastqc_out# runs fastqc on all fastq.gz n current directory. -t threads to use

# run multiqc on fastqc output folder. it automatically detects fastqc outputs
multiqc fastqc_out/* -o ./fastqc_out/multiq_fastqc --title fastqc # change title
# force interative if there are lots of files
#multiqc fastqc_out/* -o ./fastqc_out/multiqc_fastqc --interactive --title fastqc_interactive

# download and inspect output files

#==============#
# trim files using fastp
#==============#

# think about what parameters fastp may help you can change the code below

# mamba create -n fastp fastp multiqc
conda activate fastp

dir="/home/nc07/projects/metabarcoding/programs/dada2/dadaist2/data/16S_rename/"
cd $dir
ext_1="_R1.fastq.gz"
ext_2="_R2.fastq.gz"

ls $dir/*$ext_1
# for multiqc the fastp report must end in "fastp.json"

for f1 in $dir/*$ext_1
do
        f2="${f1%%$ext_1}$ext_2"
        echo $f1
        echo $f2
        fastp -i $f1 -I $f2 -o "${f1%%$ext_1}_R1.trimmed.fastq.gz" -O "${f1%%$ext_1}_R2.trimmed.fastq.gz" --max_len1 280 --max_len2 155 --report_title "${f1%%$ext_1}.report" --json "${f1%%$ext_1}.report.fastp.json" --html "${f1%%$ext_1}.report.fastp.html"
done

# run multiqc in folder to collate fastp summary outputs (does it work?)

conda activate fastqc
multiqc $dir*.json --title fastp_trimming
# multiqc *.json --title fastp_trimming_interactive --interactive

# move all trimmed files to new subdirectory
mkdir $dir/fastp_trimming # make subdirectory
mv $dir/fastp_trimming_* fastp_trimming/ # mv multiqc reports
mv $dir/*trimmed*.fastq.gz fastp_trimming/ # move trimmed reads
mv $dir/*report* fastp_trimming/ # move fastp reports

# download
