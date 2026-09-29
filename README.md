# WGS_Pf_Ethiopia_2024
Code created for analysis of Plasmodium falciparum WGS data from Ethiopia collected in 2024

## Repository overview

Scripts are listed in the order they are run in the analysis.

### 1. Data processing and variant calling

| Step | Description | Script |
|------|-------------|--------|
| 1 | Remove human reads from FASTQ files, map reads to the reference and call variants | [`fastq2vcf.sh`](fastq2vcf.sh) |
| 2 | Import per-sample VCFs into a database and merge | [`mergevcfs.sh`](mergevcfs.sh) |
| 3 | Filter the merged VCF and convert it to a genotype matrix | [`vcf2matrix.sh`](vcf2matrix.sh) |

### 2. Population genetics analyses

| Analysis | Description | Script |
|----------|-------------|--------|
| PCA | Principal component analysis of all samples | [`pca_all.qmd`](pca_all.qmd) |
| F<sub>WS</sub> | Within-host diversity (F<sub>WS</sub>) estimation | [`run_fws.sh`](run_fws.sh) |
| Ancestry | Population structure / ancestry estimation with tess3r | [`tess3r.sh`](tess3r.sh) |
| Jost's D | Pairwise population differentiation | [`process_multi_vcf_pop_jostD.qmd`](process_multi_vcf_pop_jostD.qmd) |
| Selection | Selection metrics | _TBA_ |

### 3. Statistical analysis

| Step | Description | Script |
|------|-------------|--------|
| 1 | Pre-process and filter data for regression | [`preprocess_filters.sh`](preprocess_filters.sh) |
| 2 | Regression / interaction analysis across populations (R622I) | [`multiple_pop_R622I_interaction_analysis.R`](multiple_pop_R622I_interaction_analysis.R) |
| 3 | Visualise K13 interaction results | [`K13_interaction_vis.R`](K13_interaction_vis.R) |

### 4. Map visualisations

| Description | Script |
|-------------|--------|
| Sampling location maps | [`sampling_maps.R`](sampling_maps.R) |
| Pie-chart maps of *k13*, *crt* and *hrp* markers | [`piecharts_K13_CRT_hrp.R`](piecharts_K13_CRT_hrp.R) |

