#  script to run dadaist2: https://quadram-institute-bioscience.github.io/dadaist2/

# User input: set file paths

in_dir=/home/user/path/to/read/directory/
out_dir=/home/user/path/to/read/directory/output
database=/home/nc07/path/to/database/silva_nr99_v138.1_train_set.fa.gz
meta=/home/user/path/to/metadatafile.csv # make one using dadaist2-metadata below

# example of test data
in_dir=/home/nc07/projects/metabarcoding/programs/dada2/dadaist2/data/16S_rename
out_dir=/home/nc07/projects/metabarcoding/programs/dada2/dadaist2/example-output7
database=/home/nc07/projects/metabarcoding/programs/dada2/silva_nr99_v138.1_train_set.fa.gz
meta=/home/nc07/projects/metabarcoding/programs/dada2/dadaist2/metadata2.tsv # make one using dadaist2-metadata below

# activate conda environment
conda activate dadaist2

# optionally move into project directory
cd $in_dir

# make a metadata file if one has not already been made
dadaist2-metadata -i $in_dir -o $meta

# main command - check parameters
# note - if primers not supplied switch on fastp trimming. It will skip trimming if primers are not supplied and cutadapt trimming is selected.
dadaist2 \
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
dadaist2-exporter -i $in_dir
# make a multiqc report
dadaist2-mqc-report  -i $in_dir  -o $out_dir
# find alpha diversities
dadaist2-normalize  -i $out_dir/MicrobiomeAnalyst -o OUTDIR
