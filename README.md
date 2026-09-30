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
| IBD | Identity by descent between sample pairs with hmmIBD | [`run_hmmibd_all.sh`](https://github.com/leenvh/WGS_Pf_Ethiopia_2024/blob/main/hmmibd/run_hmmibd_all.sh) |
| Tajima's D | Tajima's D across the genome | [`tajima.sh`](https://github.com/leenvh/WGS_Pf_Ethiopia_2024/blob/main/tajimasD/tajima.sh) |


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

