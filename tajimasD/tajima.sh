#!/bin/bash
# Usage:
#   cd <data_dir> && /path/to/repo/tajimasD/tajima.sh
#
# DATA_DIR: directory with metadata and where outputs are written (default: $PWD)

DATA_DIR="${DATA_DIR:-$PWD}"
cd "$DATA_DIR" || { echo "ERROR: Cannot cd to DATA_DIR=$DATA_DIR"; exit 1; }

VCF="/mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.vcf.gz"
METADATA="$DATA_DIR/Popgen_db_samples_012026_filt_FINAL_fws.csv"

# North Ethiopian districts (same definition as EHH analysis)
NORTH_DISTRICTS="Atseged Tsimbila|West Armachiho|K/Humera|Gondar|Shewa Robit Town"

# Get column indices from header
header=$(head -1 "$METADATA")
wgs_col=$(echo "$header"     | tr '\t' '\n' | grep -n "^wgs_id$"  | cut -d: -f1)
country_col=$(echo "$header" | tr '\t' '\n' | grep -n "^Country$" | cut -d: -f1)
region_col=$(echo "$header"  | tr '\t' '\n' | grep -n "^Region$"  | cut -d: -f1)
district_col=$(echo "$header"| tr '\t' '\n' | grep -n "^District$"| cut -d: -f1)
year_col=$(echo "$header"    | tr '\t' '\n' | grep -n "^Year$"     | cut -d: -f1)

echo "Column indices — wgs_id: $wgs_col | Country: $country_col | Region: $region_col | District: $district_col | Year: $year_col"
echo ""

# Helper: run vcftools TajimaD for a sample list
run_tajima() {
    local label="$1"
    local sample_file="$2"
    local outdir="$3"
    local n_samples
    n_samples=$(wc -l < "$sample_file")
    echo "  Samples: $n_samples"
    if [ "$n_samples" -eq 0 ]; then
        echo "  No samples found, skipping."
        return
    fi
    vcftools --gzvcf "$VCF" \
        --keep "$sample_file" \
        --out "${outdir}/tajimasd_${label}_10kb" \
        --TajimaD 10000
    echo "  Done: ${outdir}/tajimasd_${label}_10kb.Tajima.D"
    echo ""
}

# ------------------------------------------------------------------
# 1. By Country
# ------------------------------------------------------------------
OUTDIR="tajimasD_by_country"
mkdir -p "$OUTDIR"
countries=$(tail -n +2 "$METADATA" | awk -F'\t' -v col="$country_col" '{print $col}' | sort -u)
echo "=== By Country ==="
for group in $countries; do
    echo "Processing: $group"
    awk -F'\t' -v col_c="$country_col" -v col_s="$wgs_col" -v c="$group" \
        'NR>1 && $col_c == c {print $col_s}' "$METADATA" > "${OUTDIR}/${group}_samples.txt"
    run_tajima "$group" "${OUTDIR}/${group}_samples.txt" "$OUTDIR"
done
echo "Countries complete."
echo ""

# ------------------------------------------------------------------
# 2. By Region
# ------------------------------------------------------------------
OUTDIR="tajimasD_by_region"
mkdir -p "$OUTDIR"
regions=$(tail -n +2 "$METADATA" | awk -F'\t' -v col="$region_col" '{print $col}' | sort -u)
echo "=== By Region ==="
for group in $regions; do
    echo "Processing: $group"
    awk -F'\t' -v col_r="$region_col" -v col_s="$wgs_col" -v r="$group" \
        'NR>1 && $col_r == r {print $col_s}' "$METADATA" > "${OUTDIR}/${group}_samples.txt"
    run_tajima "$group" "${OUTDIR}/${group}_samples.txt" "$OUTDIR"
done
echo "Regions complete."
echo ""

# ------------------------------------------------------------------
# 3. Old vs New (Ethiopian samples only; Year == 2024 = "new")
# ------------------------------------------------------------------
OUTDIR="tajimasD_old_vs_new"
mkdir -p "$OUTDIR"
echo "=== Old vs New (Ethiopia only) ==="
for group in new old; do
    echo "Processing: $group"
    awk -F'\t' -v col_c="$country_col" -v col_y="$year_col" -v col_s="$wgs_col" -v g="$group" \
        'NR>1 && $col_c == "Ethiopia" && ($col_y == "2024" ? "new" : "old") == g {print $col_s}' \
        "$METADATA" > "${OUTDIR}/${group}_samples.txt"
    run_tajima "$group" "${OUTDIR}/${group}_samples.txt" "$OUTDIR"
done
echo "Old vs new complete."
echo ""

# ------------------------------------------------------------------
# 4. North vs Rest (Ethiopian samples only)
# ------------------------------------------------------------------
OUTDIR="tajimasD_north_vs_rest"
mkdir -p "$OUTDIR"
echo "=== North vs Rest (Ethiopia only) ==="
for group in north rest; do
    echo "Processing: $group"
    awk -F'\t' -v col_c="$country_col" -v col_d="$district_col" -v col_s="$wgs_col" \
        -v north="$NORTH_DISTRICTS" -v g="$group" \
        'NR>1 && $col_c == "Ethiopia" {
            is_north = ($col_d ~ north) ? "north" : "rest"
            if (is_north == g) print $col_s
        }' "$METADATA" > "${OUTDIR}/${group}_samples.txt"
    run_tajima "$group" "${OUTDIR}/${group}_samples.txt" "$OUTDIR"
done
echo "North vs rest complete."
echo ""

echo "All done."