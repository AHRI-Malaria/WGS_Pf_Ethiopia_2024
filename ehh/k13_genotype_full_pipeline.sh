#!/bin/bash
# EHH pipeline: K13 R622I genotype groups (wildtype vs homozygous mutant).
# Heterozygous (0/1) and missing samples are excluded.
#
# Usage:
#   cd <analysis_data_dir> && /path/to/repo/ehh/k13_genotype_full_pipeline.sh
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
R622I_CSV="/mnt/storage13/leen/popgen_ethiopia_WGS/WGS_Pf_Ethiopia_2024/input_files/R622I_longformat.csv"

#Extract Ethiopian samples only
bcftools view -S "$EHH_DIR/ethiopia_samples.txt" "$INPUT_VCF" -q 0.01 \
  -o 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.vqslod.filt.snps.vcf.gz -O z

#MAF 0.01
bcftools view 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.vqslod.filt.snps.vcf.gz \
  -q 0.01 -o 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.vqslod.filt.snps.vcf.gz -O z

# Calculate Tajima's D on pre-exclusion VCF (after MAF filtering), then build
# combined exclusion BED from:
#   (a) hypervarregions_032026.csv regions
#   (b) top 10% of Tajima's D windows
# Apply both exclusions together in a single bcftools step.
vcftools --gzvcf 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.vqslod.filt.snps.vcf.gz \
  --out tajimasd --TajimaD 10000

tail -n +2 tajimasd.Tajima.D | sort -k4 -rn > tajimasd_sorted.txt
n_windows=$(wc -l < tajimasd_sorted.txt)
top10pct=$(( (n_windows + 9) / 10 ))
echo "Tajima's D: $n_windows windows total, excluding top $top10pct (10%)"
head -n "$top10pct" tajimasd_sorted.txt \
  | awk '{print $1"\t"$2"\t"$2+10000}' > top10pct_tajima_windows.bed

tail -n +2 "$NEW_HYPER_VAR_CSV" | cut -f1-3 > hypervar_regions_032026.bed

cat hypervar_regions_032026.bed top10pct_tajima_windows.bed \
  | sort -k1,1 -k2,2n > combined_exclusion_032026.bed
echo "Combined exclusion BED: $(wc -l < combined_exclusion_032026.bed) regions"

echo "SNPs before exclusion:"
bcftools view -H 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.vqslod.filt.snps.vcf.gz | wc -l

bcftools view \
  --targets-file ^combined_exclusion_032026.bed \
  012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.vqslod.filt.snps.vcf.gz \
  -o 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.hypervarrem.vqslod.filt.snps.tajrem.vcf.gz \
  -O z

bcftools index 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.hypervarrem.vqslod.filt.snps.tajrem.vcf.gz

echo "SNPs after exclusion (hypervar + top 10% Tajima's D):"
bcftools view -H 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.hypervarrem.vqslod.filt.snps.tajrem.vcf.gz | wc -l

#Convert to binary matrix
python ~/software/fastq2matrix/scripts/vcf2matrix.py \
    --vcf 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.hypervarrem.vqslod.filt.snps.tajrem.vcf.gz \
    --no-iupacgt \
    --na NA \
    --threads 32

#Create annotation TSV
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ\n' \
  012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.hypervarrem.vqslod.filt.snps.tajrem.vcf.gz \
  > 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.hypervarrem.vqslod.filt.snps.tajrem.wide.tsv

awk -F'\t' 'BEGIN{OFS="\t"} {split($5,a,"|"); print $1,$2,$3,$4,a[2]}' \
  012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.hypervarrem.vqslod.filt.snps.tajrem.wide.tsv \
  > 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.hypervarrem.vqslod.filt.snps.tajrem.wide.gene.tsv

sed -i '1s/^/chr\tpos\tref\talt\tGene_name\n/' \
  012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.hypervarrem.vqslod.filt.snps.tajrem.wide.gene.tsv

#Create metadata with r622i_group column joined from R622I genotype CSV
# r622i_wt = 0/0, r622i_mut = 1/1; heterozygous (0/1) excluded
python3 << PYEOF
import pandas as pd

meta = pd.read_csv("$EHH_DIR/Popgen_db_samples_012026_filt_FINAL_fws.csv", sep="\t")
r622i = pd.read_csv("$R622I_CSV")
r622i = r622i[r622i["Genotype"].isin(["0/0", "1/1"])]
r622i["r622i_group"] = r622i["Genotype"].map({"0/0": "r622i_wt", "1/1": "r622i_mut"})
r622i_uniq = r622i[["Sample", "r622i_group"]].drop_duplicates("Sample").rename(columns={"Sample": "wgs_id"})
merged = meta.merge(r622i_uniq, on="wgs_id", how="left")
merged.to_csv("metadata_r622i_group.tsv", sep="\t", index=False)
print(merged["r622i_group"].value_counts().to_string())
PYEOF

#Run selection
conda activate selection

ANNOTATION="012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.hypervarrem.vqslod.filt.snps.tajrem.wide.gene.tsv"

printf 'r622i_wt\nr622i_mut\n' | xargs -I {} -P 2 sh -c "Rscript '$SCRIPT_DIR/prepare_input_rehh_per_category.R' \
  -d . \
  -b 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.hypervarrem.vqslod.filt.snps.tajrem.mat.bin \
  --remove_chr Pf3D7_API_v3,Pf_M76611 \
  --metadata metadata_r622i_group.tsv \
  --annotation '$ANNOTATION' \
  -c {} \
  --label_category r622i_group \
  --label_fws fws \
  --fws_th 0.95 \
  --label_id wgs_id \
  --forced_mixed"

Rscript "$SCRIPT_DIR/calculate_rehh_metrics.R" \
  -d . \
  --prefix scanned_haplotypes \
  --remove_chr Pf3D7_API_v3,Pf_M76611 \
  --list_category <(printf 'r622i_wt\nr622i_mut\n') \
  --annotation "$ANNOTATION" \
  --gene_product /mnt/storage13/leen/plasmodium_falciparum/pf_genome_product.tsv
