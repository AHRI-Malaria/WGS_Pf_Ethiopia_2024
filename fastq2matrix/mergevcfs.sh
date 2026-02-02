#!/bin/bash

######### Import samples into database

# Leucodepleted samples
python ~/software/fastq2matrix/scripts/merge_vcfs.py import \
    --sample-file /mnt/storage13/leen/leucodepleted_WGS/leucodepleted_samples.txt \
    --prefix 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen \
    --ref /mnt/storage13/leen/plasmodium_falciparum/Pfalciparum.genome.fasta \
    --vcf-dir /mnt/storage13/leen/leucodepleted_WGS/hostile \
    --threads 32

# tes
python ~/software/fastq2matrix/scripts/merge_vcfs.py import \
    --sample-file /mnt/storage13/leen/tes_WGS/tes_d0_samples.txt \
    --prefix 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen \
    --ref /mnt/storage13/leen/plasmodium_falciparum/Pfalciparum.genome.fasta \
    --vcf-dir /mnt/storage13/leen/tes_WGS/hostile \
    --threads 32


# Feleke et al
python ~/software/fastq2matrix/scripts/merge_vcfs.py import \
    --sample-file /mnt/storage13/leen/feleke_2021_WGS/RunAccList.txt \
    --prefix 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen \
    --ref /mnt/storage13/leen/plasmodium_falciparum/Pfalciparum.genome.fasta \
    --vcf-dir /mnt/storage13/leen/feleke_2021_WGS/bwa_m_flag \
    --threads 32


# mgen samples
python ~/software/fastq2matrix/scripts/merge_vcfs.py import \
    --sample-file /mnt/storage13/leen/popgen_ethiopia_WGS/072025_pf_mgen_ethiopia_popgen/mgen_database_samples.txt \
    --prefix 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen \
    --ref /mnt/storage13/leen/plasmodium_falciparum/Pfalciparum.genome.fasta \
    --vcf-dir /mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/vcfs \
    --threads 32


# Niare et al samples
python ~/software/fastq2matrix/scripts/merge_vcfs.py import \
    --sample-file /mnt/storage13/leen/niare_2025_WGS/SRR_Acc_List.txt \
    --prefix 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen \
    --ref /mnt/storage13/leen/plasmodium_falciparum/Pfalciparum.genome.fasta \
    --vcf-dir /mnt/storage13/leen/niare_2025_WGS/bwa_m_flag \
    --threads 32


# Mgen A675V SEA samples
python ~/software/fastq2matrix/scripts/merge_vcfs.py import \
    --sample-file /mnt/storage13/leen/popgen_ethiopia_WGS/112025_pf_mgen_ug_ethiopia_popgen/675_high_conf.txt \
    --prefix 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen \
    --ref /mnt/storage13/leen/plasmodium_falciparum/Pfalciparum.genome.fasta \
    --vcf-dir /mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/vcfs \
    --threads 32


# 100 most recent India Malariagen samples:
python ~/software/fastq2matrix/scripts/merge_vcfs.py import \
    --sample-file /mnt/storage13/leen/popgen_ethiopia_WGS/122025_pf_mgen_ug_india_ethiopia_popgen/india_samples_to_include.txt \
    --prefix 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen \
    --ref /mnt/storage13/leen/plasmodium_falciparum/Pfalciparum.genome.fasta \
    --vcf-dir /mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/vcfs \
    --threads 32


# Mihreteab samples:
python ~/software/fastq2matrix/scripts/merge_vcfs.py import \
    --sample-file /mnt/storage13/leen/mihreteab_2023_WGS/samples.txt \
    --prefix 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen \
    --ref /mnt/storage13/leen/plasmodium_falciparum/Pfalciparum.genome.fasta \
    --vcf-dir /mnt/storage13/leen/mihreteab_2023_WGS \
    --threads 32


######### Merge VCFs
python ~/software/fastq2matrix/scripts/merge_vcfs.py genotype --ref /mnt/storage13/leen/plasmodium_falciparum/Pfalciparum.genome.fasta --prefix 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen









