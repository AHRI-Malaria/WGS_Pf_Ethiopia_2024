#!/bin/bash
# Run hmmIBD for all comparison groupings:
#   1. By Country
#   2. By Region
#   3. Old vs New (Ethiopia only; Year 2024 = "new")
#   4. North vs Rest (Ethiopia only)
#
# Usage:
#   cd <data_dir> && /path/to/repo/hmmibd/run_hmmibd_all.sh
#
# DATA_DIR: directory with metadata/list files and where output folders are created (default: $PWD)

DATA_DIR="${DATA_DIR:-$PWD}"
cd "$DATA_DIR" || { echo "ERROR: Cannot cd to DATA_DIR=$DATA_DIR"; exit 1; }

BINARY="/mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.mat.bin"
METADATA="$DATA_DIR/Popgen_db_samples_012026_filt_FINAL_fws.csv"
RSCRIPT=~/software/malaria-hub/ibd/run_hmmIBD_per_category.R
SUMMARY_RSCRIPT=~/software/malaria-hub/ibd/summary_hmmIBD_results.R
GENE_PRODUCT="/mnt/storage13/leen/plasmodium_falciparum/pf_genome_product.tsv"
REF_INDEX="/mnt/storage13/leen/plasmodium_falciparum/Pfalciparum.genome.fasta.fai"
NORTH_DISTRICTS="Atseged Tsimbila|West Armachiho|K/Humera|Gondar|Shewa Robit Town"

# Get column indices from metadata header
header=$(head -1 "$METADATA")
wgs_col=$(echo "$header"      | tr '\t' '\n' | grep -n "^wgs_id$"   | cut -d: -f1)
country_col=$(echo "$header"  | tr '\t' '\n' | grep -n "^Country$"  | cut -d: -f1)
region_col=$(echo "$header"   | tr '\t' '\n' | grep -n "^Region$"   | cut -d: -f1)
district_col=$(echo "$header" | tr '\t' '\n' | grep -n "^District$" | cut -d: -f1)
year_col=$(echo "$header"     | tr '\t' '\n' | grep -n "^Year$"     | cut -d: -f1)

echo "Columns — wgs_id: $wgs_col | Country: $country_col | Region: $region_col | District: $district_col | Year: $year_col"
echo ""

# Helper: run summary for a completed hmmIBD output directory
run_summary() {
    local outdir="$1"
    local list_file="$2"

    echo "  Summarising: $outdir"
    Rscript "$SUMMARY_RSCRIPT" \
        -d "$outdir" \
        --list_category "$list_file" \
        --legend ibd_matrix_hap_leg.tsv \
        --gene_product "$GENE_PRODUCT" \
        --ref_index "$REF_INDEX" \
        --maf 0.01 \
        --window_size 10000 \
        --quantile_cutoff 0.95 \
        --remove_chr Pf3D7_API_v3,Pf_M76611
}

# Helper: run hmmIBD for each line in a list file, writing outputs to outdir
run_hmmibd() {
    local list_file="$1"
    local label_category="$2"
    local metadata="$3"
    local outdir="$4"

    mkdir -p "$outdir"
    cat "$list_file" | xargs -I {} -P 2 sh -c "Rscript '$RSCRIPT' \
        -d '$outdir' \
        -b '$BINARY' \
        -m '$metadata' \
        -c {} \
        --label_category '$label_category' \
        --label_fws fws \
        --fws_th 0.95 \
        --label_id wgs_id \
        --maf 0.01 \
        --na_char NA \
        --remove_chr Pf3D7_API_v3,Pf_M76611"
}

# ------------------------------------------------------------------
# 1. By Country
# ------------------------------------------------------------------
echo "=== By Country ==="
run_hmmibd "$DATA_DIR/country_list.txt" "Country" "$METADATA" "$DATA_DIR/hmmibd_by_country"
run_summary "$DATA_DIR/hmmibd_by_country" "$DATA_DIR/country_list.txt"
echo "Countries complete."
echo ""

# ------------------------------------------------------------------
# 2. By Region
# ------------------------------------------------------------------
echo "=== By Region ==="
# Generate region list dynamically from metadata
tail -n +2 "$METADATA" | awk -F'\t' -v col="$region_col" '{print $col}' | sort -u \
    > "$DATA_DIR/region_list.txt"
run_hmmibd "$DATA_DIR/region_list.txt" "Region" "$METADATA" "$DATA_DIR/hmmibd_by_region"
run_summary "$DATA_DIR/hmmibd_by_region" "$DATA_DIR/region_list.txt"
echo "Regions complete."
echo ""

# ------------------------------------------------------------------
# 3. Old vs New (Ethiopia only)
# ------------------------------------------------------------------
echo "=== Old vs New (Ethiopia only) ==="
# Create derived metadata with temporal_group column
awk -F'\t' -v col_y="$year_col" 'BEGIN{OFS="\t"} NR==1{print $0, "temporal_group"} NR>1{print $0, ($col_y=="2024" ? "new" : "old")}' \
    "$METADATA" > "$DATA_DIR/metadata_temporal_group.csv"
# Filter to Ethiopia only
awk -F'\t' -v col_c="$country_col" 'NR==1 || $col_c == "Ethiopia"' \
    "$DATA_DIR/metadata_temporal_group.csv" > "$DATA_DIR/metadata_temporal_group_ethiopia.csv"
printf "new\nold\n" > "$DATA_DIR/old_vs_new_list.txt"
run_hmmibd "$DATA_DIR/old_vs_new_list.txt" "temporal_group" \
    "$DATA_DIR/metadata_temporal_group_ethiopia.csv" "$DATA_DIR/hmmibd_old_vs_new"
run_summary "$DATA_DIR/hmmibd_old_vs_new" "$DATA_DIR/old_vs_new_list.txt"
echo "Old vs new complete."
echo ""

# ------------------------------------------------------------------
# 4. North vs Rest (Ethiopia only)
# ------------------------------------------------------------------
echo "=== North vs Rest (Ethiopia only) ==="
# Create derived metadata with region_group column
awk -F'\t' -v col_c="$country_col" -v col_d="$district_col" -v north="$NORTH_DISTRICTS" \
    'BEGIN{OFS="\t"}
     NR==1{print $0, "region_group"}
     NR>1 && $col_c=="Ethiopia" {print $0, ($col_d ~ north ? "north" : "rest")}' \
    "$METADATA" > "$DATA_DIR/metadata_region_group.csv"
printf "north\nrest\n" > "$DATA_DIR/north_vs_rest_list.txt"
run_hmmibd "$DATA_DIR/north_vs_rest_list.txt" "region_group" \
    "$DATA_DIR/metadata_region_group.csv" "$DATA_DIR/hmmibd_north_vs_rest"
run_summary "$DATA_DIR/hmmibd_north_vs_rest" "$DATA_DIR/north_vs_rest_list.txt"
echo "North vs rest complete."
echo ""

echo "All done."
