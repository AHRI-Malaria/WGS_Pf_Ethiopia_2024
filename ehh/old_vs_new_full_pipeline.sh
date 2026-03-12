#!/bin/bash
# EHH pipeline: old (pre-2024) vs new (2024) Ethiopian samples.
#
# Usage:
#   cd <analysis_data_dir> && /path/to/repo/ehh/old_vs_new_full_pipeline.sh
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

#Extract Ethiopian samples only
bcftools view -S "$EHH_DIR/ethiopia_samples.txt" "$INPUT_VCF" -q 0.01 \
  -o 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.vqslod.filt.snps.vcf.gz -O z

#MAF 0.01
bcftools view 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.vqslod.filt.snps.vcf.gz \
  -q 0.01 -o 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.vqslod.filt.snps.vcf.gz -O z

#Remove hypervariable regions
bcftools view --targets-file ^"$EHH_DIR/hyper_variable_regions.txt" \
  012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.vqslod.filt.snps.vcf.gz \
  -o 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.hypervarrem.vqslod.filt.snps.vcf.gz -O z

#Remove top 10% Tajima's D windows
vcftools --gzvcf 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.hypervarrem.vqslod.filt.snps.vcf.gz \
  --out tajimasd --TajimaD 10000

tail -n +2 tajimasd.Tajima.D | sort -k4 -rn | head -10 | awk '{print $1"\t"$2"\t"$2+10000}' > top10_tajima_windows.bed

bcftools view \
  --targets-file ^top10_tajima_windows.bed \
  012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.hypervarrem.vqslod.filt.snps.vcf.gz \
  -o 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.hypervarrem.vqslod.filt.snps.tajrem.vcf.gz \
  -O z

bcftools index 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.hypervarrem.vqslod.filt.snps.tajrem.vcf.gz

echo "SNPs before Tajima's D filtering:"
bcftools view -H 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.hypervarrem.vqslod.filt.snps.vcf.gz | wc -l
echo "SNPs after Tajima's D filtering:"
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

#Run selection
conda activate selection

#Create metadata with temporal_group column: "new" = Year 2024, "old" = all other years
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print $0, "temporal_group"} NR>1{print $0, ($20=="2024" ? "new" : "old")}' \
    "$EHH_DIR/Popgen_db_samples_012026_filt_FINAL_fws.csv" > metadata_temporal_group.csv

ANNOTATION="012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.hypervarrem.vqslod.filt.snps.tajrem.wide.gene.tsv"

cat "$EHH_DIR/old_vs_new.txt" | xargs -I {} -P 2 sh -c "Rscript '$SCRIPT_DIR/prepare_input_rehh_per_category.R' \
  -d . \
  -b 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.ethiopia.only.MAF0.01.hypervarrem.vqslod.filt.snps.tajrem.mat.bin \
  --remove_chr Pf3D7_API_v3,Pf_M76611 \
  --metadata metadata_temporal_group.csv \
  --annotation '$ANNOTATION' \
  -c {} \
  --label_category temporal_group \
  --label_fws fws \
  --fws_th 0.95 \
  --label_id wgs_id \
  --forced_mixed"

Rscript "$SCRIPT_DIR/calculate_rehh_metrics.R" \
  -d . \
  --prefix scanned_haplotypes \
  --remove_chr Pf3D7_API_v3,Pf_M76611 \
  --list_category "$EHH_DIR/old_vs_new.txt" \
  --annotation "$ANNOTATION" \
  --gene_product /mnt/storage13/leen/plasmodium_falciparum/pf_genome_product.tsv
