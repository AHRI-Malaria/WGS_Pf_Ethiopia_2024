#!/bin/bash
# EHH pipeline: country-level comparisons (India/Myanmar/Thailand excluded).
#
# Fixed version of countries_full_pipeline.sh. Differences:
#   1. No pooled MAF filter. The original ran `bcftools view -q 0.01` on all
#      countries pooled, which removed SNPs common within one country but rare
#      overall (18-68% of each country's MAF >= 0.01 SNPs), and `-q` without
#      `:minor` tests the non-reference allele rather than the minor allele.
#      MAF is now applied only per country, by prepare_input_rehh_per_category.R.
#   2. Per-country Tajima's D exclusion. The original excluded the top 10% of
#      windows from a D computed on all countries pooled (after the MAF filter),
#      which mostly reflects between-country structure. D is now computed per
#      country on its monoclonal (Fws >= FWS_TH) samples with no MAF filter, and
#      each country excludes its own top 10% windows plus the hypervariable
#      regions. Each country therefore gets its own VCF and binary matrix;
#      Rsb/XP-EHH use only SNPs present in both countries of a pair.
#
# Usage:
#   cd <analysis_data_dir> && /path/to/repo/ehh/countries_full_pipeline_fixed.sh
#
# DATA_DIR: analysis output directory (default: $PWD)
# EHH_DIR:  shared ehh data directory (default: parent of DATA_DIR)

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="${DATA_DIR:-$PWD}"
EHH_DIR="${EHH_DIR:-$(dirname "$DATA_DIR")}"

cd "$DATA_DIR" || { echo "ERROR: Cannot cd to DATA_DIR=$DATA_DIR"; exit 1; }

source /mnt/storage13/leen/miniforge3/etc/profile.d/conda.sh
conda activate fastq2matrix

INPUT_VCF="/mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.vcf.gz"
NEW_HYPER_VAR_CSV="/mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen/nsl/hypervarregions_032026.csv"
METADATA="$EHH_DIR/Popgen_db_samples_012026_filt_FINAL_fws.csv"
COUNTRY_LIST="$EHH_DIR/country_list.txt"
PREFIX="012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2"
FWS_TH=0.95

# Column indices from metadata header
header=$(head -1 "$METADATA")
wgs_col=$(echo "$header"     | tr '\t' '\n' | grep -n "^wgs_id$"  | cut -d: -f1)
country_col=$(echo "$header" | tr '\t' '\n' | grep -n "^Country$" | cut -d: -f1)
fws_col=$(echo "$header"     | tr '\t' '\n' | grep -n "^fws$"     | cut -d: -f1)
echo "Columns — wgs_id: $wgs_col | Country: $country_col | fws: $fws_col"
echo ""

#Remove India, Myanmar, Thailand (no MAF filter here — see header)
SEA_REM_VCF="${PREFIX}.sea.rem.vqslod.filt.snps.vcf.gz"
bcftools view -S ^"$EHH_DIR/India_Thailand_Myanmar.txt" "$INPUT_VCF" -o "$SEA_REM_VCF" -O z
bcftools index -f "$SEA_REM_VCF"
bcftools query -l "$SEA_REM_VCF" > vcf_samples.txt

tail -n +2 "$NEW_HYPER_VAR_CSV" | cut -f1-3 > hypervar_regions_032026.bed

#Create annotation TSV (from the pre-exclusion VCF, so it covers every country's SNPs)
ANNOTATION="${PREFIX}.sea.rem.vqslod.filt.snps.wide.gene.tsv"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ\n' "$SEA_REM_VCF" \
  > "${PREFIX}.sea.rem.vqslod.filt.snps.wide.tsv"

awk -F'\t' 'BEGIN{OFS="\t"} {split($5,a,"|"); print $1,$2,$3,$4,a[2]}' \
  "${PREFIX}.sea.rem.vqslod.filt.snps.wide.tsv" > "$ANNOTATION"

sed -i '1s/^/chr\tpos\tref\talt\tGene_name\n/' "$ANNOTATION"

# ------------------------------------------------------------------
# Per country: subset VCF -> Tajima's D -> exclusion -> binary matrix
# ------------------------------------------------------------------
while read -r country; do
  [[ -z "$country" ]] && continue
  c="${country// /_}"
  echo "=== $country ==="

  # All samples of this country in the VCF (the rehh prep script computes MAF
  # on all of them, then keeps Fws >= FWS_TH samples for the haplotypes)
  grep -Fx -f <(awk -F'\t' -v col_c="$country_col" -v col_s="$wgs_col" -v g="$country" \
      'NR>1 && $col_c == g {print $col_s}' "$METADATA") vcf_samples.txt > "${c}_samples.txt"
  # Monoclonal subset, used for Tajima's D
  grep -Fx -f <(awk -F'\t' -v col_c="$country_col" -v col_s="$wgs_col" -v col_f="$fws_col" -v th="$FWS_TH" -v g="$country" \
      'NR>1 && $col_c == g && $col_f >= th {print $col_s}' "$METADATA") vcf_samples.txt > "${c}_samples_fws.txt"
  echo "  Samples: $(wc -l < "${c}_samples.txt") (Fws >= $FWS_TH: $(wc -l < "${c}_samples_fws.txt"))"
  if [[ ! -s "${c}_samples_fws.txt" ]]; then
    echo "  No monoclonal samples, skipping."
    continue
  fi

  # Subset to the country and drop sites with no minor allele in it. This only
  # removes sites the per-country MAF filter would drop anyway. Subsetting and
  # filtering are separate commands so the count uses this country's samples.
  country_vcf="${PREFIX}.sea.rem.${c}.vqslod.filt.snps.vcf.gz"
  bcftools view -S "${c}_samples.txt" -Ou "$SEA_REM_VCF" \
    | bcftools view -c 1:minor -o "$country_vcf" -O z
  bcftools index -f "$country_vcf"

  # Tajima's D on monoclonal samples, no MAF filter
  vcftools --gzvcf "$country_vcf" --keep "${c}_samples_fws.txt" \
    --out "tajimasd_${c}" --TajimaD 10000
  if [[ ! -s "tajimasd_${c}.Tajima.D" ]]; then
    echo "  ERROR: vcftools produced no Tajima's D for $country, skipping."
    continue
  fi

  tail -n +2 "tajimasd_${c}.Tajima.D" | awk '$4 != "nan" && $4 != "-nan"' | sort -k4,4gr \
    > "tajimasd_${c}_sorted.txt"
  n_windows=$(wc -l < "tajimasd_${c}_sorted.txt")
  top10pct=$(( (n_windows + 9) / 10 ))
  echo "  Tajima's D: $n_windows windows total, excluding top $top10pct (10%)"
  head -n "$top10pct" "tajimasd_${c}_sorted.txt" \
    | awk '{print $1"\t"$2"\t"$2+10000}' > "top10pct_tajima_windows_${c}.bed"

  cat hypervar_regions_032026.bed "top10pct_tajima_windows_${c}.bed" \
    | sort -k1,1 -k2,2n > "combined_exclusion_032026_${c}.bed"
  echo "  Combined exclusion BED: $(wc -l < "combined_exclusion_032026_${c}.bed") regions"

  final_vcf="${PREFIX}.sea.rem.${c}.hypervarrem.vqslod.filt.snps.tajrem.vcf.gz"
  echo "  SNPs before exclusion: $(bcftools view -H "$country_vcf" | wc -l)"
  bcftools view \
    --targets-file ^"combined_exclusion_032026_${c}.bed" \
    "$country_vcf" \
    -o "$final_vcf" \
    -O z
  bcftools index -f "$final_vcf"
  echo "  SNPs after exclusion (hypervar + top 10% Tajima's D): $(bcftools view -H "$final_vcf" | wc -l)"

  #Convert to binary matrix
  python ~/software/fastq2matrix/scripts/vcf2matrix.py \
      --vcf "$final_vcf" \
      --no-iupacgt \
      --na NA \
      --threads 32
  echo ""
done < "$COUNTRY_LIST"

#Run selection
conda activate selection

prep_country() {
  local country="$1"
  local c="${country// /_}"
  local matrix="${PREFIX}.sea.rem.${c}.hypervarrem.vqslod.filt.snps.tajrem.mat.bin"
  if [[ ! -s "$matrix" ]]; then
    echo "ERROR: no binary matrix for $country ($matrix), skipping rehh input."
    return 1
  fi
  Rscript "$SCRIPT_DIR/prepare_input_rehh_per_category.R" \
    -d . \
    -b "$matrix" \
    --remove_chr Pf3D7_API_v3,Pf_M76611 \
    --metadata "$METADATA" \
    --annotation "$ANNOTATION" \
    -c "$country" \
    --label_category Country \
    --label_fws fws \
    --fws_th "$FWS_TH" \
    --label_id wgs_id \
    --maf 0.01 \
    --forced_mixed
}
export -f prep_country
export SCRIPT_DIR PREFIX METADATA ANNOTATION FWS_TH

grep -v '^$' "$COUNTRY_LIST" | xargs -P 2 -I {} bash -c 'prep_country "$1"' _ {}

Rscript "$SCRIPT_DIR/calculate_rehh_metrics.R" \
  -d . \
  --prefix scanned_haplotypes \
  --remove_chr Pf3D7_API_v3,Pf_M76611 \
  --list_category "$COUNTRY_LIST" \
  --annotation "$ANNOTATION" \
  --gene_product /mnt/storage13/leen/plasmodium_falciparum/pf_genome_product.tsv
