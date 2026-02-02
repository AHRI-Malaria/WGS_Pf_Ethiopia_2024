#!/bin/bash

# Filtering VCF
python ~/software/fastq2matrix/scripts/filter_merged_vcf.py --merged-file 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.2026_01_27.genotyped.vcf.gz --prefix 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt --ref /mnt/storage13/leen/plasmodium_falciparum/Pfalciparum.genome.fasta --bqsr-vcf /mnt/storage13/leen/plasmodium_falciparum/3d7_hb3.combined.final.vcf.gz,/mnt/storage13/leen/plasmodium_falciparum/7g8_gb4.combined.final.vcf.gz,/mnt/storage13/leen/plasmodium_falciparum/hb3_dd2.combined.final.vcf.gz --include-region /mnt/storage13/leen/plasmodium_falciparum/Core_genome_Pf3D7_v3_ext.bed --vqslod 0 --missing-sample-cutoff 0.2 --cutoff-mix-GT 0.8 --gff-file /mnt/storage13/leen/plasmodium_falciparum/Pfalciparum.genome.modified.new.gff3 --threads 32 

# conversion to matrix

 