#!/bin/bash
folder=$PWD"/run_example/"
query=$folder"query.vcf"
ground_truth=$folder"ground_truth.vcf"
bed=$folder"high_conf.bed"
bam=$folder"input.bam"
fasta=$folder"chr_22.fa"
genome_bed=$folder"chr_22.fa.fai"
annovar=""
out=$folder"outdir/"
caller="TVC"

file_id="1ZDlbBrsZhfGqfLEbQIl1WRqfBC0Eh1QI"  # Replace with the Google Drive file ID
bam="run_example/input.bam"  # Replace with the desired output file name



echo "Downloading input.bam and input.bam.bai files..."




# Download the file
curl -c cookies.txt -s -L "https://drive.google.com/uc?export=download&id=${file_id}" > /dev/null
curl -Lb cookies.txt "https://drive.google.com/uc?export=download&confirm=$(awk '/download/ {print $NF}' cookies.txt)&id=${file_id}" -o "${bam}"
# Cleanup
rm cookies.txt

echo "File downloaded: ${bam}"

file_id="1GhQOGVNQLy5C0zdDRgstVoFSb4EkdGhs"  # Replace with the Google Drive file ID
bai="run_example/input.bam.bai"  # Replace with the desired output file name

# Download the file
curl -c cookies.txt -s -L "https://drive.google.com/uc?export=download&id=${file_id}" > /dev/null
curl -Lb cookies.txt "https://drive.google.com/uc?export=download&confirm=$(awk '/download/ {print $NF}' cookies.txt)&id=${file_id}" -o "${bai}"


# Cleanup
rm cookies.txt

echo "File downloaded: ${bai}"


# Check if annovar is empty
if [ -z "$annovar" ]; then
  echo "Please provide the annovar folder to RecallME v.0.2"
  exit 1  # Exit the script with an error code
fi



echo "Running RecallME.py on example files..."


python $PWD/RecallME.py \
-q $query \
-g $ground_truth \
--high_conf_bed $bed \
--query_format VCF \
--gt_format VCF \
-f $fasta \
--genome_bed $genome_bed \
-b $bam \
-a $annovar \
-o $out \
--caller $caller
