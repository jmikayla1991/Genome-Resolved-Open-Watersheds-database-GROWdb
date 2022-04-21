#!/bin/bash

#$1 list of samples corresponding to directory names

for element in $(<$1)
do
cd "$element"
mkdir raw_reads
mkdir trimmed_reads
mkdir assembly
cd Raw_Data

gunzip -c "$element"_interleaved.fastq.gz > "$element"_interleaved.fastq

/ORG-Data/scripts/deinterleave_fastq.sh < "$element"_interleaved.fastq "$element"_R1.fastq "$element"_R2.fastq

sickle pe -f "$element"_R1.fastq  -r "$element"_R2.fastq   -t sanger -o "$element"_R1_trimmed.fastq -p "$element"_R2_trimmed.fastq -s discared_R1R2.fastq

mv *trimmed.fastq ../trimmed_reads

cd ../assembly

mkdir megahit

cd megahit

megahit -1 ../../trimmed_reads/"$element"_R1_trimmed.fastq -2 ../../trimmed_reads/"$element"_R2_trimmed.fastq --k-min 31 --k-max 121 --k-step 10 -m 0.3  -t 25

cd ../../../

done
echo "reads have been successfuly trimmed"
