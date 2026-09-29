#!/bin/bash

#Necessary filtering steps prior to running fws:
#1. bi-allelic SNPs only
#2. coding regions only. Use input_files/Pf3D7_R.coding.regions.bed
#1. maf filter

#Split metadata to populations
library(readr)
library(dplyr)

met <- readr::read_csv('metadata.csv') # Load metadata
suffix <- "pf_meta" # Add additional suffix after phenotype name

met %>%
  group_by(Country) %>%
  dplyr::select(wgs_id, Country) %>%
  group_walk(~ write.table(.x, paste0(.y$Country, suffix, ".txt"), quote = F, row.names = F, col.names = F))


#Split genotyped merged VCF
ls *pf_meta.txt > population_files_list.txt


cat population_files_list.txt | xargs -I {} sh -c 'bcftools view -c1 -S {} -Oz --threads 20 -o {}.genotyped.vcf.gz myvcf.vcf.gz'

#Run moimix for each population
ls *.genotyped.vcf.gz | sed 's/.genotyped.vcf.gz//g' > population_names_list.txt


cat population_names_list.txt | xargs -I {} sh -c 'Rscript calculate_fws.R -d . -f {}.genotyped.vcf.gz'

