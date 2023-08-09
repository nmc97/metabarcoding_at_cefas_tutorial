# run cutadapt on multiple samples
# untested
# Auther NC

# Set primer sequences
FWD=GTGCCAGCAGYYGCGGTAATAC
REV=CACRACACGAGCTGACGACA
FWD.RC=GTATTACCGCRRCTGCTGGCAC
REV.RC=TGTCGTCAGCTCGTGTYGTG

# Set paths and extensions
dir="/path/to/directory/" # set path to files
cd $dir # move into that directory
ext_1="_R1.fastq.gz" # set file extension for forward reads
ext_2="_R2.fastq.gz" # set file extension for reverse reads

ls $dir/*$ext_1 # check files exist with the variables set above

# run a loop
# loop will run through each sample and print each line of code need to trim using cutadapt
# check that these lines of code look correct. then remove `echo` and the quotes around the line of code and run the loop again.
for f1 in $dir/*$ext_1
do
        f2="${f1%%$ext_1}$ext_2" # this uses the list of forward reads to build the file names of reverse reads
        echo $f1 # reads out the forward read file name
        echo $f2 # reads out the reverse read file name
        # if you wish to check the loop is finding the correct files then write  `ls $f1` and `ls $f2` and comment out the fastp command.
        fastp -i $f1 -I $f2 -o "${f1%%$ext_1}_R1.trimmed.fastq.gz" -O "${f1%%$ext_1}_R2.trimmed.fastq.gz" --max_len1 280 --max_len2 155 --report_title "${f1%%$ext_1}.report" --json "${f1%%$ext_1}.report.fastp.json" --html "${f1%%$ext_1}.report.fastp.html"
	echo "cutadapt -g FWD -a REV.RC -G REV -A FWD.RC -n 2 --only-trimmed \
-o orientation_1/${f1%%$ext_1}_forward.amplicon${ext1} -p orientation_1/${f1%%$ext_1}_forward.amplicon${ext2} \
${f1%%$ext_1}${ext_1} ${f1%%$ext_1}${ext2}"
done 