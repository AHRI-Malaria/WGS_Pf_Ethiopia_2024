#!/bin/bash

# Constants
FASTQ_DIR="dir"
REF="path/Pfalciparum.genome.fasta"
BQSR_VCFS="path/3d7_hb3.combined.final.vcf.gz,path/7g8_gb4.combined.final.vcf.gz,path/hb3_dd2.combined.final.vcf.gz"

for R1 in "${FASTQ_DIR}"/*WG24*_1.clean_1.fastq.gz; do
    # Determine R2 filename
    R2="${R1/_1.clean_1.fastq.gz/_2.clean_2.fastq.gz}"

    # Extract prefix (e.g., SRR15000240)
    PREFIX=$(basename "$R1" | sed -n 's/.*\(WG24[^_]*\)_lib.*/\1/p')

    # Check if R2 file exists
    if [ ! -f "$R2" ]; then
        echo "Warning: Paired file for $R1 not found. Skipping."
        continue
    fi

    echo "Processing sample: $PREFIX"
    python /mnt/storage13/leen/software/fastq2matrix/scripts/fastq2vcf.py all \
        --read1 "$R1" \
        --read2 "$R2" \
        --ref "$REF" \
        --prefix "$PREFIX" \
        --bqsr-vcf "$BQSR_VCFS" \
        --threads 32
done