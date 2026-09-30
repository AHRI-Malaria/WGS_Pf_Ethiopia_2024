# WGS_Pf_Ethiopia_2024
Code created for analysis of Plasmodium falciparum WGS data from Ethiopia collected in 2024

## Repository overview

Scripts are listed in the order they are run in the analysis.

### 1. Data processing and variant calling

| Step | Description | Script |
|------|-------------|--------|
| 1 | Remove human reads from FASTQ files, map reads to the reference and call variants | [`fastq2vcf.sh`](https://github.com/leenvh/WGS_Pf_Ethiopia_2024/blob/main/fastq2matrix/fastq2vcf.sh) |
| 2 | Import per-sample VCFs into a database and merge | [`mergevcfs.sh`](https://github.com/leenvh/WGS_Pf_Ethiopia_2024/blob/main/fastq2matrix/mergevcfs.sh) |
| 3 | Filter the merged VCF and convert it to a genotype matrix | [`vcf2matrix.sh`](https://github.com/leenvh/WGS_Pf_Ethiopia_2024/blob/main/fastq2matrix/vcf2matrix.sh) |

### 2. Population genetics analyses

| Analysis | Description | Script |
|----------|-------------|--------|
| PCA | Principal component analysis of all samples | [`pca_all.qmd`](https://github.com/leenvh/WGS_Pf_Ethiopia_2024/blob/main/pca/pca_all.qmd) |
| F<sub>WS</sub> | Within-host diversity (F<sub>WS</sub>) estimation | [`run_fws.sh`](https://github.com/leenvh/WGS_Pf_Ethiopia_2024/blob/main/fws/run_fws.sh) |
| Ancestry | Population structure / ancestry estimation with tess3r | [`tess3r.sh`](https://github.com/leenvh/WGS_Pf_Ethiopia_2024/blob/main/tess3r/tess3r.sh) |
| Jost's D | Pairwise population differentiation | [`process_multi_vcf_pop_jostD.qmd`](https://github.com/leenvh/WGS_Pf_Ethiopia_2024/blob/main/josts%20D/process_multi_vcf_pop_jostD.qmd) |

### 3. Selection metrics

| Analysis | Description | Script |
|----------|-------------|--------|
| EHH | Extended haplotype homozygosity, one script per population comparison | [`ehh/`](https://github.com/leenvh/WGS_Pf_Ethiopia_2024/tree/main/ehh) |
| EHH decay | EHH decay plots for SNPs of interest | [`ehh/ehh_plots/`](https://github.com/leenvh/WGS_Pf_Ethiopia_2024/tree/main/ehh/ehh_plots) |
| IBD | Identity by descent between sample pairs with hmmIBD for specified comparison groups| [`run_hmmibd_all.sh`](https://github.com/leenvh/WGS_Pf_Ethiopia_2024/blob/main/hmmibd/run_hmmibd_all.sh) |
| Tajima's D | Tajima's D across the genome for specified comparison groups | [`tajima.sh`](https://github.com/leenvh/WGS_Pf_Ethiopia_2024/blob/main/tajimasD/tajima.sh) |
| Composite selection signals | Multi-tiered analytical framework combining all four selection metrics for specified comparison groups | [`comp_selection/`](https://github.com/leenvh/WGS_Pf_Ethiopia_2024/tree/main/comp_selection) |


### 4. Statistical analysis

| Step | Description | Script |
|------|-------------|--------|
| 1 | Pre-process and filter data for regression | [`preprocess_filters.sh`](https://github.com/leenvh/WGS_Pf_Ethiopia_2024/blob/main/regression/preprocess_filters.sh) |
| 2 | Regression / interaction analysis across populations (R622I) | [`multiple_pop_R622I_interaction_analysis.R`](https://github.com/leenvh/WGS_Pf_Ethiopia_2024/blob/main/regression/multiple_pop_R622I_interaction_analysis.R) |
| 3 | Visualise K13 interaction results | [`K13_interaction_vis.R`](https://github.com/leenvh/WGS_Pf_Ethiopia_2024/blob/main/regression/K13_interaction_vis.R) |

### 5. Map visualisations

| Description | Script |
|-------------|--------|
| Sampling location maps | [`sampling_maps.R`](https://github.com/leenvh/WGS_Pf_Ethiopia_2024/blob/main/maps/sampling_maps.R) |
| Pie-chart maps of *k13*, *crt* and *hrp* markers | [`piecharts_K13_CRT_hrp.R`](https://github.com/leenvh/WGS_Pf_Ethiopia_2024/blob/main/maps/piecharts_K13_CRT_hrp.R) |


## Requirements
### Bash environment
Bash scripts scripts rely on functioning installations of bwa, python, samtools, gatk, and bcftools.<br />
The analysis was done with the following versions: bwa version 0.7.18-r1243, python version 3.12.11, samtools version 1.22, gatk v4.1.4.1, bcftools version 1.22 <br />


### R environment
R scripts were executed in R version 4.4.2 using the following packages: ape (5.8.1), pegas (1.3), msa (1.38.0), dplyr (1.1.4), psych (2.4.5), officer (0.6.7), flextable (0.9.7), tidyr (1.3.1), readxl (1.4.5), purrr (1.0.4), writexl (1.5.1), reshape2 (1.4.4), tidyverse (2.0.0), ggplot2 (3.5.2), scales (1.4.0), sf (1.0.20), grid (4.4.2), stringr (1.5.1), ggforce (0.4.2), grid (4.4.2), tmap (4.0), cowplot (1.1.3), scatterpie (0.2.5), gridExtra (2.3), ggrepel (0.9.6), data.table (1.18.2.1), rehh (3.2.3) , vcfR (1.15.0), qqman (0.1.9), R.utils (1.8.2)

### Python environment
Python scripts were run in python version 3.12.7 and using the following packages: numpy (1.26.4), scikit-allel (1.3.13), pysam (0.23.3), plotly (5.24.1), scipy (1.14.1), matplotlib (3.9.3), pandas (2.2.3), tqdm (4.66.5)
