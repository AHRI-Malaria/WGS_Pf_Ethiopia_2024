#!/bin/bash
# EHH pipeline: pfcrt K76 vs 76T within the co-endemic sites
# (Gondar, Mizan, Metehera Fentale; Country == Ethiopia & Year == 2024).
#
# GROUP DEFINITION
# The split is made on pfcrt K76T alone (Pf3D7_07_v3:403,625), taken from
# input_files/K76T_longformat.csv:
#     K76 = 0/0 (wild-type)      76T = 1/1 (mutant)
# Heterozygous (0/1), missing (./.) and second-alternate calls are excluded.
# Codons 75 and 76 are multiallelic and therefore absent from the .bi. VCF used
# downstream, which is why the genotype comes from the long-format CSV.
#
# The groups are named for the codon-76 genotype rather than the CVIET/CVMNK
# haplotype because the other haplotype-defining codons cannot be resolved
# here: M74I (403,620) has cohort AC=6 and is 0/0 in every co-endemic-site
# sample, and no variant is called anywhere in codon 75 (403,621-403,623).
#
# Group sizes at fws >= 0.95: 76T 32, K76 23 (of 55 monoclonal / 73 total).
#
#
# Usage:
#   cd <analysis_data_dir> && /path/to/repo/ehh/pfcrt_genotype_full_pipeline_fixed.sh
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
K76T_CSV="/mnt/storage13/leen/popgen_ethiopia_WGS/WGS_Pf_Ethiopia_2024/input_files/K76T_longformat.csv"
METADATA="$EHH_DIR/Popgen_db_samples_012026_filt_FINAL_fws.csv"
FWS_TH=0.95

# Districts making up the co-endemic set (as written in the metadata)
CO_ENDEMIC="Gondar|Mizan|Metehera Fentale"

PREFIX="012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.pfcrt.coendemic.only"

# Get column indices from the metadata header (tab-separated despite the .csv name)
header=$(head -1 "$METADATA")
wgs_col=$(echo "$header"      | tr '\t' '\n' | grep -n "^wgs_id$"   | cut -d: -f1)
country_col=$(echo "$header"  | tr '\t' '\n' | grep -n "^Country$"  | cut -d: -f1)
district_col=$(echo "$header" | tr '\t' '\n' | grep -n "^District$" | cut -d: -f1)
year_col=$(echo "$header"     | tr '\t' '\n' | grep -n "^Year$"     | cut -d: -f1)

echo "Columns — wgs_id: $wgs_col | Country: $country_col | District: $district_col | Year: $year_col"
echo ""

#Build the co-endemic-site sample list, restricted to samples present in the VCF
awk -F'\t' -v col_s="$wgs_col" -v col_c="$country_col" -v col_y="$year_col" -v col_d="$district_col" \
    -v co="$CO_ENDEMIC" \
    'NR>1 && $col_c=="Ethiopia" && $col_y=="2024" && $col_d ~ co {print $col_s}' \
    "$METADATA" | sort -u > pfcrt_samples_all.txt
bcftools query -l "$INPUT_VCF" | sort -u > vcf_samples.txt
comm -12 pfcrt_samples_all.txt vcf_samples.txt > pfcrt_samples.txt
echo "Co-endemic-site samples in metadata: $(wc -l < pfcrt_samples_all.txt)"
echo "Co-endemic-site samples present in VCF: $(wc -l < pfcrt_samples.txt)"
if [ ! -s pfcrt_samples.txt ]; then
  echo "ERROR: no co-endemic-site samples found in VCF"; exit 1
fi

#Extract the co-endemic-site samples only (no MAF filter — see header)
BASE_VCF="$PREFIX.vqslod.filt.snps.vcf.gz"
bcftools view -S pfcrt_samples.txt "$INPUT_VCF" -o "$BASE_VCF" -O z
bcftools index -f "$BASE_VCF"
bcftools query -l "$BASE_VCF" | sort -u > base_vcf_samples.txt

tail -n +2 "$NEW_HYPER_VAR_CSV" | cut -f1-3 > hypervar_regions_032026.bed

#Create annotation TSV (from the pre-exclusion VCF, so it covers every group's SNPs)
ANNOTATION="$PREFIX.vqslod.filt.snps.wide.gene.tsv"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ\n' "$BASE_VCF" \
  > "$PREFIX.vqslod.filt.snps.wide.tsv"

awk -F'\t' 'BEGIN{OFS="\t"} {split($5,a,"|"); print $1,$2,$3,$4,a[2]}' \
  "$PREFIX.vqslod.filt.snps.wide.tsv" > "$ANNOTATION"

sed -i '1s/^/chr\tpos\tref\talt\tGene_name\n/' "$ANNOTATION"

#Create metadata with a pfcrt_group column joined from the K76T genotype CSV,
#restricted to the co-endemic sites. Genotypes are phased ("1|1") for part of
#the file, so separators are normalised before filtering.
python3 << PYEOF
import pandas as pd

meta = pd.read_csv("$METADATA", sep="\t")
co = ["Gondar", "Mizan", "Metehera Fentale"]
in_co = (meta["Country"] == "Ethiopia") & (meta["Year"].astype(str) == "2024") \
        & (meta["District"].isin(co))

k76t = pd.read_csv("$K76T_CSV")
k76t["gt"] = k76t["Genotype"].str.replace("|", "/", regex=False)
k76t = k76t[k76t["gt"].isin(["0/0", "1/1"])]
k76t["pfcrt_group"] = k76t["gt"].map({"0/0": "K76", "1/1": "76T"})
k76t_uniq = k76t[["Sample", "pfcrt_group"]].drop_duplicates("Sample").rename(
    columns={"Sample": "wgs_id"})

merged = meta.merge(k76t_uniq, on="wgs_id", how="left")
# Only samples from the co-endemic sites keep a group label
merged.loc[~in_co.reindex(merged.index, fill_value=False), "pfcrt_group"] = pd.NA
merged.to_csv("metadata_pfcrt_group.tsv", sep="\t", index=False)

print(merged["pfcrt_group"].value_counts().to_string())
print("monoclonal (fws >= 0.95):")
print(merged[merged["fws"] >= 0.95]["pfcrt_group"].value_counts().to_string())
PYEOF

GROUP_METADATA="metadata_pfcrt_group.tsv"
LABEL="pfcrt_group"
GROUP_LIST="pfcrt_group_list.txt"
printf 'K76\n76T\n' > "$GROUP_LIST"

# ------------------------------------------------------------------
# Per group: subset VCF -> Tajima's D -> exclusion -> binary matrix
# ------------------------------------------------------------------
# Prints wgs_id for rows whose LABEL column equals grp; with fws_only=1 only
# rows with Fws >= th (a missing Fws never passes)
GROUP_SAMPLES_AWK='NR==1{for(i=1;i<=NF;i++){if($i=="wgs_id")cs=i; if($i==lab)cl=i; if($i=="fws")cf=i}; next}
  $cl == grp && (!fws_only || ($cf != "" && $cf != "NA" && $cf+0 >= th)) {print $cs}'

# The list is read on fd 3 so commands inside the loop cannot consume it
while read -r group <&3; do
  [[ -z "$group" ]] && continue
  g="${group// /_}"
  echo "=== $group ==="

  # All samples of this group in the base VCF (the rehh prep script computes MAF
  # on all of them, then keeps Fws >= FWS_TH samples for the haplotypes), and
  # the monoclonal subset used for Tajima's D
  awk -F'\t' -v lab="$LABEL" -v grp="$group" -v th="$FWS_TH" -v fws_only=0 "$GROUP_SAMPLES_AWK" \
    "$GROUP_METADATA" | sort -u | comm -12 - base_vcf_samples.txt > "${g}_samples.txt"
  awk -F'\t' -v lab="$LABEL" -v grp="$group" -v th="$FWS_TH" -v fws_only=1 "$GROUP_SAMPLES_AWK" \
    "$GROUP_METADATA" | sort -u | comm -12 - base_vcf_samples.txt > "${g}_samples_fws.txt"
  echo "  Samples: $(wc -l < "${g}_samples.txt") (Fws >= $FWS_TH: $(wc -l < "${g}_samples_fws.txt"))"
  if [[ ! -s "${g}_samples_fws.txt" ]]; then
    echo "  No monoclonal samples, skipping."
    continue
  fi

  # Subset to the group and drop sites with no minor allele in it. This only
  # removes sites the per-group MAF filter would drop anyway. Subsetting and
  # filtering are separate commands so the count uses this group's samples.
  group_vcf="${PREFIX}.${g}.vqslod.filt.snps.vcf.gz"
  bcftools view -S "${g}_samples.txt" -Ou "$BASE_VCF" \
    | bcftools view -c 1:minor -o "$group_vcf" -O z
  bcftools index -f "$group_vcf"

  # Tajima's D on monoclonal samples, no MAF filter
  vcftools --gzvcf "$group_vcf" --keep "${g}_samples_fws.txt" \
    --out "tajimasd_${g}" --TajimaD 10000
  if [[ ! -s "tajimasd_${g}.Tajima.D" ]]; then
    echo "  ERROR: vcftools produced no Tajima's D for $group, skipping."
    continue
  fi

  tail -n +2 "tajimasd_${g}.Tajima.D" | awk '$4 != "nan" && $4 != "-nan"' | sort -k4,4gr \
    > "tajimasd_${g}_sorted.txt"
  n_windows=$(wc -l < "tajimasd_${g}_sorted.txt")
  top10pct=$(( (n_windows + 9) / 10 ))
  echo "  Tajima's D: $n_windows windows total, excluding top $top10pct (10%)"
  head -n "$top10pct" "tajimasd_${g}_sorted.txt" \
    | awk '{print $1"\t"$2"\t"$2+10000}' > "top10pct_tajima_windows_${g}.bed"

  cat hypervar_regions_032026.bed "top10pct_tajima_windows_${g}.bed" \
    | sort -k1,1 -k2,2n > "combined_exclusion_032026_${g}.bed"
  echo "  Combined exclusion BED: $(wc -l < "combined_exclusion_032026_${g}.bed") regions"

  final_vcf="${PREFIX}.${g}.hypervarrem.vqslod.filt.snps.tajrem.vcf.gz"
  echo "  SNPs before exclusion: $(bcftools view -H "$group_vcf" | wc -l)"
  bcftools view \
    --targets-file ^"combined_exclusion_032026_${g}.bed" \
    "$group_vcf" \
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
done 3< "$GROUP_LIST"

#Run selection
conda activate selection

prep_group() {
  local group="$1"
  local g="${group// /_}"
  local matrix="${PREFIX}.${g}.hypervarrem.vqslod.filt.snps.tajrem.mat.bin"
  if [[ ! -s "$matrix" ]]; then
    echo "ERROR: no binary matrix for $group ($matrix), skipping rehh input."
    return 1
  fi
  Rscript "$SCRIPT_DIR/prepare_input_rehh_per_category.R" \
    -d . \
    -b "$matrix" \
    --remove_chr Pf3D7_API_v3,Pf_M76611 \
    --metadata "$GROUP_METADATA" \
    --annotation "$ANNOTATION" \
    -c "$group" \
    --label_category "$LABEL" \
    --label_fws fws \
    --fws_th "$FWS_TH" \
    --label_id wgs_id \
    --maf 0.01 \
    --forced_mixed
}
export -f prep_group
export SCRIPT_DIR PREFIX GROUP_METADATA LABEL ANNOTATION FWS_TH

grep -v '^$' "$GROUP_LIST" | xargs -P 2 -I {} bash -c 'prep_group "$1"' _ {}

Rscript "$SCRIPT_DIR/calculate_rehh_metrics.R" \
  -d . \
  --prefix scanned_haplotypes \
  --remove_chr Pf3D7_API_v3,Pf_M76611 \
  --list_category "$GROUP_LIST" \
  --annotation "$ANNOTATION" \
  --gene_product /mnt/storage13/leen/plasmodium_falciparum/pf_genome_product.tsv
