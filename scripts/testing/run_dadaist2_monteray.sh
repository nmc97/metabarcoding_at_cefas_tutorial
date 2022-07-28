#  script to run dadaist2: https://quadram-institute-bioscience.github.io/dadaist2/

# User input: set file paths
in_dir=/home/user/path/to/read/directory/
out_dir=/home/user/path/to/read/directory/output
database=/home/nc07/path/to/database/silva_nr99_v138.1_train_set.fa.gz
meta=/home/user/path/to/metadatafile.csv # make one using dadaist2-metadata below

# example of test data
in_dir=/home/nc07/projects/metabarcoding/workshop/monterey/16S/fastq
out_dir=/home/nc07/projects/metabarcoding/workshop/monterey/outputs/dadaist2
database=/home/nc07/projects/metabarcoding/programs/dada2/silva_nr99_v138.1_train_set.fa.gz
meta=/home/nc07/projects/metabarcoding/workshop/monterey/16S/fastq/metadata2.tsv # make one using dadaist2-metadata below

# activate conda environment
conda activate dadaist2

# optionally move into project directory
cd $in_dir

# make a metadata file if one has not already been made
dadaist2-metadata -i $in_dir -o $meta --for-tag "_1" --rev-tag "_2"

# main command - check parameters
# note - if primers not supplied switch on fastp trimming. It will skip trimming if primers are not supplied and cutadapt trimming is selected.
dadaist2 \
-input-directory $in_dir  \
-output-directory $out_dir \
-database $database \
-metadata $meta \
-threads 12 \
-trim-primer-for 10 \
-trim-primer-rev 10 \
-trunc-len-1 150 \
-trunc-len-2 150 \
-min-qual 28 \
-maxee1 2 \
-maxee2 2 \
-save-rds \
-fastp \
-1 _1 \
-2 _2 \
-verbose

# export to get MetagenomeAnalyist compatable files
dadaist2-exporter -i $out_dir
# make a multiqc report
dadaist2-mqc-report  -i $out_dir  -o $out_dir/multiqc
# find alpha diversities
dadaist2-normalize  -i $out_dir/feature-table.tsv -o $out_dir/normalise
