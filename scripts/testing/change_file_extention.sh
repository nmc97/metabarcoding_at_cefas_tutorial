
dir=$1
ext=$2

files=$dir/*$ext
for file in $files; do
  echo $file

done
