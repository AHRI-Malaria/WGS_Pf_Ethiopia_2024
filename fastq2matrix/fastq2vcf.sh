#!/bin/bash

# Cleaning newly sequenced samples with hostile

# leucodepleted samples
cat ../leucodepleted_samples.txt | xargs -I {} -P 10 sh -c \
'R1=$(find ../fastq -name "*{}_lib*_1.fastq.gz"); R2=$(find ../fastq -name "*{}_lib*_2.fastq.gz"); [ -n "$R1" ] && [ -n "$R2" ] && hostile clean --fastq1 $R1 --fastq2 $R2 --threads 8 --force || echo "Files not found for {}"'

# TES samples
ls *_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//' > samples_R1R2.txt 
ls *_1.fastq.gz | sed 's/_1.fastq.gz//' > samples_12.txt

cat ../fastq/samples_R1R2.txt | xargs -I {} -P 10 sh -c \ 'hostile clean --fastq1 ../fastq/{}_R1_001.fastq.gz --fastq2 ../fastq/{}_R2_001.fastq.gz --threads 8 --force' 
cat ../fastq/samples_12.txt | xargs -I {} -P 10 sh -c \ 'hostile clean --fastq1 ../fastq/{}_1.fastq.gz --fastq2 ../fastq/{}_2.fastq.gz --threads 8'


# Mapping and variant calling
./run_fastq2vcf.sh

