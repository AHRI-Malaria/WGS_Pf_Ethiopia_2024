#!/usr/bin/env bash
# =============================================================================
# PRE-PROCESSING FOR EPISTASIS ANALYSIS
# Run this script before multiple_pop_R622I_interaction_analysis.R
#
# Scenarios processed:
#   EthEri   — Ethiopia + Eritrea (steps 1-7, already done)
#   HOA      — Ethiopia + Eritrea + Sudan (PCA only, VCF already done)
#   Ethiopia — Ethiopia only
#   EthNew   — Ethiopia 2024 only
#
# Steps per scenario:
#   1. Filter for target samples
#   2. Remove monomorphic sites
#   3. Filter to monoclonal samples (Fws >= 0.95)
#   4. Apply MAF >= 0.01 filter
#   5. Convert to matrix
#   6. Make annotation file
#   7. Extract PCA components
# =============================================================================

set -euo pipefail

WORKDIR="/mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen/k13"
FWS_META="/mnt/storage13/leen/popgen_ethiopia_WGS/WGS_Pf_Ethiopia_2024/input_files/Popgen_db_samples_012026_filt_FINAL_fws.csv"
VCF_IN="${WORKDIR}/../012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.vcf.gz"
VCF_HOA="${WORKDIR}/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.EthEri.vcf.gz"
VCF_FILT="${WORKDIR}/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.EthEri.filt.vcf.gz"
VCF_MONO="${WORKDIR}/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.EthEri.filt.mono.vcf.gz"
VCF_OUT="${WORKDIR}/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.EthEri.filt.mono.maf0.01.vcf.gz"
MONO_SAMPLES="${WORKDIR}/monoclonal_samples.txt"
PREFIX="012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps"
THREADS=16

cd "${WORKDIR}"

# =============================================================================
# STEP 1 — Filter for Ethiopia and Eritrea samples only
# =============================================================================
bcftools view -S Eth_Eri_samples.txt \
  "${VCF_IN}" \
  -Oz -o "${VCF_HOA}"

bcftools index "${VCF_HOA}"

# =============================================================================
# STEP 2 — Remove monomorphic sites
# =============================================================================
echo "Removing monomorphic sites..."

bcftools view \
    --min-ac 1:minor \
    --output-type z \
    --threads "${THREADS}" \
    --output "${VCF_FILT}" \
    "${VCF_HOA}"

bcftools index "${VCF_FILT}"

echo "  SNPs before monomorphic filter: $(bcftools view -H "${VCF_HOA}"  | wc -l)"
echo "  SNPs after  monomorphic filter: $(bcftools view -H "${VCF_FILT}" | wc -l)"

# =============================================================================
# STEP 3 — Extract monoclonal sample IDs (Fws >= 0.95)
# =============================================================================
echo "Extracting monoclonal samples (Fws >= 0.95)..."

# wgs_id is column 1, fws is the last column (NF)
awk -F',' 'NR > 1 && $NF >= 0.95 { print $1 }' "${FWS_META}" > "${MONO_SAMPLES}"

echo "  Monoclonal samples found: $(wc -l < "${MONO_SAMPLES}")"

echo "Subsetting VCF to monoclonal samples..."

bcftools view \
    --samples-file "${MONO_SAMPLES}" \
    --output-type z \
    --force-samples \
    --threads "${THREADS}" \
    --output "${VCF_MONO}" \
    "${VCF_FILT}"

bcftools index "${VCF_MONO}"

echo "  Samples in output: $(bcftools query -l "${VCF_MONO}" | wc -l)"

# =============================================================================
# STEP 4 — Apply MAF >= 0.01 filter
# =============================================================================
echo "Applying MAF >= 0.01 filter..."

bcftools view \
    --min-af 0.01:minor \
    --output-type z \
    --threads "${THREADS}" \
    --output "${VCF_OUT}" \
    "${VCF_MONO}"

bcftools index "${VCF_OUT}"

echo "  SNPs before MAF filter: $(bcftools view -H "${VCF_MONO}" | wc -l)"
echo "  SNPs after  MAF filter: $(bcftools view -H "${VCF_OUT}"  | wc -l)"

# =============================================================================
# STEP 5 — Convert to matrix for R script
# =============================================================================
echo "Converting to genotype matrix..."

python ~/software/fastq2matrix/scripts/vcf2matrix.py \
    --vcf "${VCF_OUT}" \
    --no-iupacgt \
    --na NA \
    --threads "${THREADS}"

echo ""
echo "Done. Update nina_interaction_analysis.R to use:"
echo "  geno_file <- \"$(basename "${VCF_OUT}" .vcf.gz).mat.bin\""


# =============================================================================
# STEP 6 — Extract and parse annotations
# =============================================================================
echo "Extracting annotations from VCF..."

bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%BCSQ\n' \
    "${VCF_OUT}" > annotations.tsv

# Parse BCSQ field (pipe-delimited): field1|gene_name1|gene_name2|...
# Output: chr, pos, ref, alt, Gene_name (formatted as "gene_name1 (gene_name2)")
awk 'BEGIN { OFS="\t"; print "chr\tpos\tref\talt\tGene_name" } {
    n = split($5, a, "|")
    gene = (n >= 3) ? a[2] " (" a[3] ")" : $5
    print $1, $2, $3, $4, gene
}' annotations.tsv > annotations_final.tsv

echo "  Annotations written to annotations_final.tsv"

# =============================================================================
# STEP 7 — Compute PCA and format coordinates
# =============================================================================
echo "Running PCA..."

plink \
    --vcf "${VCF_OUT}" \
    --pca 10 \
    --double-id \
    --allow-extra-chr \
    --out PC_032026

# Convert eigenvec to CSV (columns: FID IID PC1..PC10 → keep IID as sample)
awk 'BEGIN { OFS=","; print "sample,PC1,PC2,PC3,PC4,PC5" } \
     { print $2, $3, $4, $5, $6, $7 }' \
    PC_032026.eigenvec > PC_032026_pca_coordinates.csv

echo "  PCA coordinates written to PC_032026_pca_coordinates.csv"

# =============================================================================
# STEP 8 — Generate HOA PCA (VCF already exists)
# =============================================================================
echo "Running PCA for HOA..."

HOA_VCF_OUT="${WORKDIR}/${PREFIX}.HOA.filt.mono.maf0.01.vcf.gz"

plink \
    --vcf "${HOA_VCF_OUT}" \
    --pca 10 \
    --double-id \
    --allow-extra-chr \
    --out PC_HOA

awk 'BEGIN { OFS=","; print "sample,PC1,PC2,PC3,PC4,PC5" } \
     { print $2, $3, $4, $5, $6, $7 }' \
    PC_HOA.eigenvec > PC_HOA_pca_coordinates.csv

echo "  PCA coordinates written to PC_HOA_pca_coordinates.csv"

# =============================================================================
# STEP 9 — Generate sample lists for Ethiopia and EthNew
# =============================================================================
echo "Generating sample lists..."

awk -F',' 'NR>1 && $13=="Ethiopia" { print $1 }' \
    "${FWS_META}" > "${WORKDIR}/Ethiopia_samples.txt"
echo "  Ethiopia samples: $(wc -l < "${WORKDIR}/Ethiopia_samples.txt")"

awk -F',' 'NR>1 && $13=="Ethiopia" && $20=="2024" { print $1 }' \
    "${FWS_META}" > "${WORKDIR}/EthNew_samples.txt"
echo "  EthNew (2024) samples: $(wc -l < "${WORKDIR}/EthNew_samples.txt")"

# =============================================================================
# STEP 10 — Run pipeline for Ethiopia and EthNew scenarios
# =============================================================================

run_scenario() {
    local SCENARIO="$1"
    local SAMPLES_FILE="${WORKDIR}/${SCENARIO}_samples.txt"

    local VCF_S="${WORKDIR}/${PREFIX}.${SCENARIO}.vcf.gz"
    local VCF_S_FILT="${WORKDIR}/${PREFIX}.${SCENARIO}.filt.vcf.gz"
    local VCF_S_MONO="${WORKDIR}/${PREFIX}.${SCENARIO}.filt.mono.vcf.gz"
    local VCF_S_OUT="${WORKDIR}/${PREFIX}.${SCENARIO}.filt.mono.maf0.01.vcf.gz"

    echo ""
    echo "=== Processing scenario: ${SCENARIO} ==="

    # Step 1: subset to target samples
    echo "  Subsetting samples..."
    bcftools view -S "${SAMPLES_FILE}" \
        "${VCF_IN}" \
        -Oz -o "${VCF_S}"
    bcftools index "${VCF_S}"

    # Step 2: remove monomorphic sites
    echo "  Removing monomorphic sites..."
    bcftools view \
        --min-ac 1:minor \
        --output-type z \
        --threads "${THREADS}" \
        --output "${VCF_S_FILT}" \
        "${VCF_S}"
    bcftools index "${VCF_S_FILT}"
    echo "    SNPs before: $(bcftools view -H "${VCF_S}"      | wc -l)"
    echo "    SNPs after:  $(bcftools view -H "${VCF_S_FILT}" | wc -l)"

    # Step 3: filter to monoclonal samples (Fws >= 0.95)
    echo "  Filtering to monoclonal samples..."
    bcftools view \
        --samples-file "${MONO_SAMPLES}" \
        --output-type z \
        --force-samples \
        --threads "${THREADS}" \
        --output "${VCF_S_MONO}" \
        "${VCF_S_FILT}"
    bcftools index "${VCF_S_MONO}"
    echo "    Samples retained: $(bcftools query -l "${VCF_S_MONO}" | wc -l)"

    # Step 4: MAF >= 0.01 filter
    echo "  Applying MAF >= 0.01..."
    bcftools view \
        --min-af 0.01:minor \
        --output-type z \
        --threads "${THREADS}" \
        --output "${VCF_S_OUT}" \
        "${VCF_S_MONO}"
    bcftools index "${VCF_S_OUT}"
    echo "    SNPs before: $(bcftools view -H "${VCF_S_MONO}" | wc -l)"
    echo "    SNPs after:  $(bcftools view -H "${VCF_S_OUT}"  | wc -l)"

    # Step 5: convert to matrix
    echo "  Converting to genotype matrix..."
    python ~/software/fastq2matrix/scripts/vcf2matrix.py \
        --vcf "${VCF_S_OUT}" \
        --no-iupacgt \
        --na NA \
        --threads "${THREADS}"

    # Step 6: PCA
    echo "  Running PCA..."
    plink \
        --vcf "${VCF_S_OUT}" \
        --pca 10 \
        --double-id \
        --allow-extra-chr \
        --out "PC_${SCENARIO}"

    awk 'BEGIN { OFS=","; print "sample,PC1,PC2,PC3,PC4,PC5" } \
         { print $2, $3, $4, $5, $6, $7 }' \
        "PC_${SCENARIO}.eigenvec" > "PC_${SCENARIO}_pca_coordinates.csv"

    echo "  Done. Matrix: $(basename "${VCF_S_OUT}" .vcf.gz).mat.bin"
    echo "  PCA:          PC_${SCENARIO}_pca_coordinates.csv"
}

run_scenario "Ethiopia"
run_scenario "EthNew"

echo ""
echo "All scenarios complete."
