#!/bin/bash
# Run hmmIBD for all comparison groupings:
#   1. By Country
#   2. By Region
#   3. Old vs New (Ethiopia only; Year 2024 = "new")
#   4. North vs Rest (Ethiopia only)
#   5. Kelch13 R622I: WT (0/0) vs mutant (1/1) — heterozygous (0/1) excluded
#   6. HRP status: wild vs hrp3_only vs mutant (hrp2-/hrp3-)  — hrp2- only excluded
#   7. By Site (Ethiopia 2024 collection sites: Gondar, Gambella, Metehera
#      Fentale, Mizan, Arba Minch)
#   8. Sites grouping (Co_endemic_sites = Gondar + Mizan + Metehera Fentale,
#      vs Gambella vs Arba Minch)
#   9. pfcrt K76 vs 76T within the co-endemic sites
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
R622I_CSV="/mnt/storage13/leen/popgen_ethiopia_WGS/WGS_Pf_Ethiopia_2024/input_files/R622I_longformat.csv"
HRP_META="/mnt/storage13/leen/popgen_ethiopia_WGS/WGS_Pf_Ethiopia_2024/input_files/sample_meta_with_hrp_status_updated.csv"
K76T_CSV="/mnt/storage13/leen/popgen_ethiopia_WGS/WGS_Pf_Ethiopia_2024/input_files/K76T_longformat.csv"
# Co-endemic 2024 collection sites (same definition as the EHH and Tajima's D pipelines)
CO_ENDEMIC="Gondar|Mizan|Metehera Fentale"
PF_SIZES="/mnt/storage13/leen/popgen_ethiopia_WGS/072025_pf_mgen_ethiopia_popgen/hmmIBD/Ethiopia/pf.sizes"

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
        --legend "$outdir/ibd_matrix_hap_leg.tsv" \
        --gene_product "$GENE_PRODUCT" \
        --ref_index "$REF_INDEX" \
        --maf 0.01 \
        --window_size 10000 \
        --quantile_cutoff 0.95 \
        --remove_chr Pf3D7_API_v3,Pf_M76611
}

# Centromere BED lines appended to every chromosome painting BED.
# chr10 is an estimate; all others from PMC5052046.
CENTROMERE_LINES="1\t457932\t460311\tcentromere
2\t447301\t450450\tcentromere
3\t597817\t600275\tcentromere
4\t642004\t644529\tcentromere
5\t455741\t457252\tcentromere
6\t478653\t480971\tcentromere
7\t809246\t811716\tcentromere
8\t299080\t301403\tcentromere
9\t1242138\t1244483\tcentromere
10\t1100000\t1102200\tcentromere
11\t831969\t834245\tcentromere
12\t1282774\t1285067\tcentromere
13\t1168128\t1170425\tcentromere
14\t1071524\t1075089\tcentromere"

# Helper: chromosome painting from a completed hmmIBD summary directory.
#   Produces one painting per category found in the *_q0.95.tsv (col 7).
paint_chromosomes() {
    local outdir="$1"
    local title_prefix="$2"

    local qtsv
    qtsv=$(ls "$outdir"/*_q0.95.tsv 2>/dev/null | head -1)
    if [[ -z "$qtsv" ]]; then
        echo "  WARNING: no *_q0.95.tsv found in $outdir — skipping chromosome painting."
        return
    fi

    # Get unique categories from col 7
    local grps
    grps=$(awk 'BEGIN{FS="\t"} NR>1 {print $7}' "$qtsv" | sort -u)

    for grp in $grps; do
        echo "  Chromosome painting: $title_prefix — $grp"
        local bedfile="$outdir/ibd_chrom_paint_${grp}.bed"

        # Filter to this category, convert to BED (0-based start)
        awk -v g="$grp" 'BEGIN{FS=OFS="\t"} NR>1 && $7==g {print $1, $2-1, $3, "high_IBD"}' \
            "$qtsv" > "$bedfile"

        # Append centromere positions
        printf "%b\n" "$CENTROMERE_LINES" >> "$bedfile"

        # Run jcvi (jcvi lives in the fastq2matrix conda env)
        ( cd "$outdir" && /mnt/storage13/leen/miniforge3/envs/fastq2matrix/bin/python \
            -m jcvi.graphics.chromosome "ibd_chrom_paint_${grp}.bed" \
            --size="$PF_SIZES" \
            --title="$title_prefix — $grp" \
            --gauge )
    done
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
paint_chromosomes "$DATA_DIR/hmmibd_by_country" "High IBD regions — by country"
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
paint_chromosomes "$DATA_DIR/hmmibd_by_region" "High IBD regions — by region"
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
paint_chromosomes "$DATA_DIR/hmmibd_old_vs_new" "High IBD regions — old vs new (Ethiopia)"
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
paint_chromosomes "$DATA_DIR/hmmibd_north_vs_rest" "High IBD regions — north vs rest (Ethiopia)"
echo "North vs rest complete."
echo ""

# ------------------------------------------------------------------
# 5. Kelch13 R622I: WT (0/0) vs mutant (1/1)
#    Heterozygous (0/1) and missing (./.) samples are excluded.
#    R622I_longformat.csv columns: Chrom,Position,Ref,Alt,Sample,Genotype
# ------------------------------------------------------------------
echo "=== Kelch13 R622I: WT vs Mutant ==="
# Build Sample -> kelch13_r622i lookup (exclude 0/1 and ./.)
awk -F',' 'NR>1 && ($6=="0/0" || $6=="1/1") {print $5 "\t" ($6=="0/0" ? "wt" : "mut")}' \
    "$R622I_CSV" > "$DATA_DIR/r622i_lookup.tmp"
# Join with main metadata: restrict to Ethiopia and only keep rows whose wgs_id appears in the lookup
awk -F'\t' -v col_id="$wgs_col" -v col_c="$country_col" \
    'BEGIN{OFS="\t"; while((getline line < ARGV[1])>0){split(line,a,"\t"); grp[a[1]]=a[2]} ARGV[1]=""}
     NR==1{print $0, "kelch13_r622i"}
     NR>1 && $col_c=="Ethiopia" && ($col_id in grp){print $0, grp[$col_id]}' \
    "$DATA_DIR/r622i_lookup.tmp" "$METADATA" > "$DATA_DIR/metadata_r622i.csv"
rm "$DATA_DIR/r622i_lookup.tmp"
printf "wt\nmut\n" > "$DATA_DIR/r622i_list.txt"
run_hmmibd "$DATA_DIR/r622i_list.txt" "kelch13_r622i" \
    "$DATA_DIR/metadata_r622i.csv" "$DATA_DIR/hmmibd_r622i"
run_summary "$DATA_DIR/hmmibd_r622i" "$DATA_DIR/r622i_list.txt"
paint_chromosomes "$DATA_DIR/hmmibd_r622i" "High IBD regions — Kelch13 R622I (WT vs mutant)"
echo "R622I complete."
echo ""

# ------------------------------------------------------------------
# 6. HRP status: wild vs hrp3_only vs mutant (hrp2-/hrp3-)
#    hrp2_only isolates are excluded from this comparison.
#    HRP metadata columns: Sample,...,hrp_status (col 8)
# ------------------------------------------------------------------
echo "=== HRP status: wild vs hrp3_only vs mutant ==="
# Build Sample -> hrp_status lookup (keep only the 3 groups of interest)
awk -F',' 'NR>1 && ($8=="wild" || $8=="hrp3_only" || $8=="mutant") {print $1 "\t" $8}' \
    "$HRP_META" > "$DATA_DIR/hrp_lookup.tmp"
# Join with main metadata: restrict to Ethiopia and only keep rows whose wgs_id appears in the lookup
awk -F'\t' -v col_id="$wgs_col" -v col_c="$country_col" \
    'BEGIN{OFS="\t"; while((getline line < ARGV[1])>0){split(line,a,"\t"); grp[a[1]]=a[2]} ARGV[1]=""}
     NR==1{print $0, "hrp_status"}
     NR>1 && $col_c=="Ethiopia" && ($col_id in grp){print $0, grp[$col_id]}' \
    "$DATA_DIR/hrp_lookup.tmp" "$METADATA" > "$DATA_DIR/metadata_hrp.csv"
rm "$DATA_DIR/hrp_lookup.tmp"
printf "wild\nhrp3_only\nmutant\n" > "$DATA_DIR/hrp_list.txt"
run_hmmibd "$DATA_DIR/hrp_list.txt" "hrp_status" \
    "$DATA_DIR/metadata_hrp.csv" "$DATA_DIR/hmmibd_hrp_status"
run_summary "$DATA_DIR/hmmibd_hrp_status" "$DATA_DIR/hrp_list.txt"
paint_chromosomes "$DATA_DIR/hmmibd_hrp_status" "High IBD regions — HRP status (wild vs hrp3_only vs mutant)"
echo "HRP status complete."
echo ""

# ------------------------------------------------------------------
# 7. By Site (Ethiopian 2024 collection sites)
#    Gondar, Gambella, Metehera Fentale, Mizan, Arba Minch.
#    Site labels are underscore-normalised (e.g. "Arba Minch" -> "Arba_Minch")
#    so they survive shell word-splitting and filename construction.
# ------------------------------------------------------------------
echo "=== By Site (Ethiopia 2024) ==="
# Create derived metadata with a 'site' column, restricted to the 2024 Ethiopian samples
awk -F'\t' -v col_c="$country_col" -v col_y="$year_col" -v col_d="$district_col" \
    'BEGIN{OFS="\t"}
     NR==1{print $0, "site"}
     NR>1 && $col_c=="Ethiopia" && $col_y=="2024" && $col_d != "" {
        s=$col_d; gsub(/ +/, "_", s); print $0, s
     }' "$METADATA" > "$DATA_DIR/metadata_site.csv"
# Generate site list dynamically from the derived metadata (last column)
awk -F'\t' 'NR>1 {print $NF}' "$DATA_DIR/metadata_site.csv" | sort -u \
    > "$DATA_DIR/site_list.txt"
echo "  Sites: $(tr '\n' ' ' < "$DATA_DIR/site_list.txt")"
run_hmmibd "$DATA_DIR/site_list.txt" "site" \
    "$DATA_DIR/metadata_site.csv" "$DATA_DIR/hmmibd_by_site"
run_summary "$DATA_DIR/hmmibd_by_site" "$DATA_DIR/site_list.txt"
paint_chromosomes "$DATA_DIR/hmmibd_by_site" "High IBD regions — by site (Ethiopia 2024)"
echo "Sites complete."
echo ""

# ------------------------------------------------------------------
# 8. Sites grouping (Ethiopia 2024)
#    Co_endemic_sites = Gondar + Mizan + Metehera Fentale, compared with
#    Gambella and Arba Minch. Gambella and Arba_Minch are recomputed here
#    rather than reused from section 7 so the comparison is self-contained;
#    downstream scripts tell the two apart on the 'comparison' column.
# ------------------------------------------------------------------
echo "=== Sites grouping (Ethiopia 2024) ==="
# Create derived metadata with a 'site_group' column, restricted to the 2024
# Ethiopian samples
awk -F'\t' -v col_c="$country_col" -v col_y="$year_col" -v col_d="$district_col" \
    -v co="$CO_ENDEMIC" \
    'BEGIN{OFS="\t"}
     NR==1{print $0, "site_group"}
     NR>1 && $col_c=="Ethiopia" && $col_y=="2024" && $col_d != "" {
        if ($col_d ~ co) { print $0, "Co_endemic_sites" }
        else { s=$col_d; gsub(/ +/, "_", s); print $0, s }
     }' "$METADATA" > "$DATA_DIR/metadata_site_group.csv"
awk -F'\t' 'NR>1 {print $NF}' "$DATA_DIR/metadata_site_group.csv" | sort -u \
    > "$DATA_DIR/sites_grouping_list.txt"
echo "  Groups: $(tr '\n' ' ' < "$DATA_DIR/sites_grouping_list.txt")"
run_hmmibd "$DATA_DIR/sites_grouping_list.txt" "site_group" \
    "$DATA_DIR/metadata_site_group.csv" "$DATA_DIR/hmmibd_sites_grouping"
run_summary "$DATA_DIR/hmmibd_sites_grouping" "$DATA_DIR/sites_grouping_list.txt"
paint_chromosomes "$DATA_DIR/hmmibd_sites_grouping" "High IBD regions — sites grouping (Ethiopia 2024)"
echo "Sites grouping complete."
echo ""

# ------------------------------------------------------------------
# 9. pfcrt K76 vs 76T (co-endemic sites only)
#    Groups are defined by K76T (Pf3D7_07_v3:403,625):
#      K76 = 0/0 (wild-type), 76T = 1/1 (mutant)
#    Heterozygous, missing and second-alternate calls are excluded. K76T is
#    multiallelic and absent from the .bi. VCF, so genotypes come from
#    K76T_longformat.csv; part of that file is phased ("1|1"), so separators
#    are normalised first, and its quoted Alt field ("C,G") means the CSV needs
#    a quote-aware parser rather than -F','.
#    The groups are named for the codon-76 genotype rather than the CVIET/CVMNK
#    haplotype because M74I and N75E are not variant in this callset — see
#    ehh/pfcrt_genotype_full_pipeline.sh for the detail.
# ------------------------------------------------------------------
echo "=== pfcrt: K76 vs 76T (co-endemic sites) ==="
python3 -c '
import csv, sys
with open(sys.argv[1], newline="") as fh:
    for row in csv.DictReader(fh):
        gt = row["Genotype"].replace("|", "/")
        if gt == "0/0":
            print(row["Sample"] + "\tK76")
        elif gt == "1/1":
            print(row["Sample"] + "\t76T")
' "$K76T_CSV" > "$DATA_DIR/pfcrt_lookup.tmp"
# Join with main metadata: co-endemic sites only, and only rows with a call
awk -F'\t' -v col_id="$wgs_col" -v col_c="$country_col" -v col_y="$year_col" \
    -v col_d="$district_col" -v co="$CO_ENDEMIC" \
    'BEGIN{OFS="\t"; while((getline line < ARGV[1])>0){split(line,a,"\t"); grp[a[1]]=a[2]} ARGV[1]=""}
     NR==1{print $0, "pfcrt_group"}
     NR>1 && $col_c=="Ethiopia" && $col_y=="2024" && $col_d ~ co && ($col_id in grp){print $0, grp[$col_id]}' \
    "$DATA_DIR/pfcrt_lookup.tmp" "$METADATA" > "$DATA_DIR/metadata_pfcrt.csv"
rm "$DATA_DIR/pfcrt_lookup.tmp"
printf "K76\n76T\n" > "$DATA_DIR/pfcrt_list.txt"
run_hmmibd "$DATA_DIR/pfcrt_list.txt" "pfcrt_group" \
    "$DATA_DIR/metadata_pfcrt.csv" "$DATA_DIR/hmmibd_pfcrt"
run_summary "$DATA_DIR/hmmibd_pfcrt" "$DATA_DIR/pfcrt_list.txt"
paint_chromosomes "$DATA_DIR/hmmibd_pfcrt" "High IBD regions — pfcrt K76 vs 76T (co-endemic sites)"
echo "pfcrt complete."
echo ""

# ------------------------------------------------------------------
# 10. Combine win10kb_annotated.tsv files from all comparisons
# ------------------------------------------------------------------
echo "=== Combining win10kb_annotated.tsv files ==="
combined="$DATA_DIR/hmmibd_combined_win10kb_annotated.tsv"
header_written=0
> "$combined"

for f in "$DATA_DIR"/hmmibd_*/*win10kb_annotated.tsv; do
    [[ -f "$f" ]] || continue
    comparison=$(basename "$(dirname "$f")" | sed 's/^hmmibd_//')
    if [[ $header_written -eq 0 ]]; then
        head -1 "$f" | awk '{print $0 "\tcomparison"}' >> "$combined"
        header_written=1
    fi
    tail -n +2 "$f" | awk -v cmp="$comparison" '{print $0 "\t" cmp}' >> "$combined"
done

if [[ $header_written -eq 1 ]]; then
    echo "  Combined file: $combined"
    wc -l "$combined"
else
    echo "  WARNING: no win10kb_annotated.tsv files found — nothing combined."
fi
echo ""

echo "All done."
