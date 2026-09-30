#!/bin/bash
# Usage:
#   cd <data_dir> && /path/to/repo/tajimasD/tajima.sh
#
# DATA_DIR: directory with metadata and where outputs are written (default: $PWD)

DATA_DIR="${DATA_DIR:-$PWD}"
cd "$DATA_DIR" || { echo "ERROR: Cannot cd to DATA_DIR=$DATA_DIR"; exit 1; }

VCF="/mnt/storage13/leen/popgen_ethiopia_WGS/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen/012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.vcf.gz"
METADATA="$DATA_DIR/Popgen_db_samples_012026_filt_FINAL_fws.csv"
R622I_CSV="/mnt/storage13/leen/popgen_ethiopia_WGS/WGS_Pf_Ethiopia_2024/input_files/R622I_longformat.csv"
HRP_META="/mnt/storage13/leen/popgen_ethiopia_WGS/WGS_Pf_Ethiopia_2024/input_files/sample_meta_with_hrp_status_updated.csv"
K76T_CSV="/mnt/storage13/leen/popgen_ethiopia_WGS/WGS_Pf_Ethiopia_2024/input_files/K76T_longformat.csv"
FWS_TH=0.95

# vcftools lives in the fastq2matrix conda env. Reference it by explicit path so
# the script works from a plain shell (e.g. a fresh screen session) without the
# env activated first — otherwise every vcftools call fails with "command not
# found", no .Tajima.D files are written, and the script carries on to annotate
# and merge using whatever stale output happens to be on disk.
VCFTOOLS="${VCFTOOLS:-/mnt/storage13/leen/miniforge3/envs/fastq2matrix/bin/vcftools}"
if [ ! -x "$VCFTOOLS" ]; then
    VCFTOOLS=$(command -v vcftools)
fi
if [ -z "$VCFTOOLS" ]; then
    echo "ERROR: vcftools not found. Set VCFTOOLS=/path/to/vcftools or activate the fastq2matrix env."
    exit 1
fi
echo "Using vcftools: $VCFTOOLS"

# North Ethiopian districts (same definition as EHH analysis)
NORTH_DISTRICTS="Atseged Tsimbila|West Armachiho|K/Humera|Gondar|Shewa Robit Town"

# Co-endemic 2024 collection sites, pooled into one group in section 8 and
# forming the sample set of section 9 (same definition as the EHH pipelines)
CO_ENDEMIC="Gondar|Mizan|Metehera Fentale"

# Get column indices from header
header=$(head -1 "$METADATA")
wgs_col=$(echo "$header"     | tr '\t' '\n' | grep -n "^wgs_id$"  | cut -d: -f1)
country_col=$(echo "$header" | tr '\t' '\n' | grep -n "^Country$" | cut -d: -f1)
region_col=$(echo "$header"  | tr '\t' '\n' | grep -n "^Region$"  | cut -d: -f1)
district_col=$(echo "$header"| tr '\t' '\n' | grep -n "^District$"| cut -d: -f1)
year_col=$(echo "$header"    | tr '\t' '\n' | grep -n "^Year$"     | cut -d: -f1)
fws_col=$(echo "$header"     | tr '\t' '\n' | grep -n "^fws$"      | cut -d: -f1)

echo "Column indices — wgs_id: $wgs_col | Country: $country_col | Region: $region_col | District: $district_col | Year: $year_col | fws: $fws_col"
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
    "$VCFTOOLS" --gzvcf "$VCF" \
        --keep "$sample_file" \
        --maf 0.01 \
        --out "${outdir}/tajimasd_${label}_10kb" \
        --TajimaD 10000
    if [ ! -s "${outdir}/tajimasd_${label}_10kb.Tajima.D" ]; then
        echo "  ERROR: vcftools produced no output for $label"
        return 1
    fi
    echo "  Done: ${outdir}/tajimasd_${label}_10kb.Tajima.D"
    echo ""
}

# ------------------------------------------------------------------
# 1. By Country
# ------------------------------------------------------------------
OUTDIR="tajimasD_by_country"
mkdir -p "$OUTDIR"
countries=$(tail -n +2 "$METADATA" | awk -F'\t' -v col="$country_col" -v f="$fws_col" -v th="$FWS_TH" \
    'NR>1 && $f >= th {print $col}' | sort -u)
echo "=== By Country ==="
for group in $countries; do
    echo "Processing: $group"
    awk -F'\t' -v col_c="$country_col" -v col_s="$wgs_col" -v col_f="$fws_col" -v th="$FWS_TH" -v c="$group" \
        'NR>1 && $col_c == c && $col_f >= th {print $col_s}' "$METADATA" > "${OUTDIR}/${group}_samples.txt"
    run_tajima "$group" "${OUTDIR}/${group}_samples.txt" "$OUTDIR"
done
echo "Countries complete."
echo ""

# ------------------------------------------------------------------
# 2. By Region
# ------------------------------------------------------------------
OUTDIR="tajimasD_by_region"
mkdir -p "$OUTDIR"
regions=$(tail -n +2 "$METADATA" | awk -F'\t' -v col="$region_col" -v f="$fws_col" -v th="$FWS_TH" \
    'NR>1 && $f >= th {print $col}' | sort -u)
echo "=== By Region ==="
for group in $regions; do
    echo "Processing: $group"
    awk -F'\t' -v col_r="$region_col" -v col_s="$wgs_col" -v col_f="$fws_col" -v th="$FWS_TH" -v r="$group" \
        'NR>1 && $col_r == r && $col_f >= th {print $col_s}' "$METADATA" > "${OUTDIR}/${group}_samples.txt"
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
    awk -F'\t' -v col_c="$country_col" -v col_y="$year_col" -v col_s="$wgs_col" \
        -v col_f="$fws_col" -v th="$FWS_TH" -v g="$group" \
        'NR>1 && $col_c == "Ethiopia" && ($col_y == "2024" ? "new" : "old") == g && $col_f >= th {print $col_s}' \
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
        -v col_f="$fws_col" -v th="$FWS_TH" -v north="$NORTH_DISTRICTS" -v g="$group" \
        'NR>1 && $col_c == "Ethiopia" && $col_f >= th {
            is_north = ($col_d ~ north) ? "north" : "rest"
            if (is_north == g) print $col_s
        }' "$METADATA" > "${OUTDIR}/${group}_samples.txt"
    run_tajima "$group" "${OUTDIR}/${group}_samples.txt" "$OUTDIR"
done
echo "North vs rest complete."
echo ""

# ------------------------------------------------------------------
# 5. Kelch13 R622I: WT (0/0) vs mutant (1/1)
#    Heterozygous (0/1) and missing (./.) samples are excluded.
#    Only monoclonal samples (FWS >= FWS_TH) retained.
#    R622I_longformat.csv columns: Chrom,Position,Ref,Alt,Sample,Genotype
# ------------------------------------------------------------------
OUTDIR="tajimasD_r622i"
mkdir -p "$OUTDIR"
echo "=== Kelch13 R622I: WT vs Mutant ==="
awk -F',' 'NR>1 && ($6=="0/0" || $6=="1/1") {print $5 "\t" ($6=="0/0" ? "wt" : "mut")}' \
    "$R622I_CSV" > "$OUTDIR/r622i_lookup.tmp"
for group in wt mut; do
    echo "Processing: $group"
    awk -F'\t' -v col_s="$wgs_col" -v col_c="$country_col" -v col_f="$fws_col" -v th="$FWS_TH" -v g="$group" \
        'BEGIN{while((getline line < ARGV[1])>0){split(line,a,"\t"); grp[a[1]]=a[2]} ARGV[1]=""}
         NR>1 && $col_c=="Ethiopia" && ($col_s in grp) && grp[$col_s]==g && $col_f >= th {print $col_s}' \
        "$OUTDIR/r622i_lookup.tmp" "$METADATA" > "${OUTDIR}/${group}_samples.txt"
    run_tajima "$group" "${OUTDIR}/${group}_samples.txt" "$OUTDIR"
done
rm "$OUTDIR/r622i_lookup.tmp"
echo "R622I complete."
echo ""

# ------------------------------------------------------------------
# 6. HRP status: wild vs hrp3_only vs mutant (hrp2-/hrp3-)
#    hrp2_only isolates are excluded from this comparison.
#    Only monoclonal samples (FWS >= FWS_TH) retained.
#    HRP metadata columns: Sample,...,hrp_status (col 8)
#    Note: gsub(/\r/,"") handles Windows CRLF line endings in the HRP file.
# ------------------------------------------------------------------
OUTDIR="tajimasD_hrp_status"
mkdir -p "$OUTDIR"
echo "=== HRP status: wild vs hrp3_only vs mutant ==="
awk -F',' 'NR>1 {gsub(/\r/, "", $8); if ($8=="wild" || $8=="hrp3_only" || $8=="mutant") print $1 "\t" $8}' \
    "$HRP_META" > "$OUTDIR/hrp_lookup.tmp"
for group in wild hrp3_only mutant; do
    echo "Processing: $group"
    awk -F'\t' -v col_s="$wgs_col" -v col_c="$country_col" -v col_f="$fws_col" -v th="$FWS_TH" -v g="$group" \
        'BEGIN{while((getline line < ARGV[1])>0){split(line,a,"\t"); grp[a[1]]=a[2]} ARGV[1]=""}
         NR>1 && $col_c=="Ethiopia" && ($col_s in grp) && grp[$col_s]==g && $col_f >= th {print $col_s}' \
        "$OUTDIR/hrp_lookup.tmp" "$METADATA" > "${OUTDIR}/${group}_samples.txt"
    run_tajima "$group" "${OUTDIR}/${group}_samples.txt" "$OUTDIR"
done
rm "$OUTDIR/hrp_lookup.tmp"
echo "HRP status complete."
echo ""

# ------------------------------------------------------------------
# 7. By Site (Ethiopian 2024 collection sites)
#    Gondar, Gambella, Metehera Fentale, Mizan, Arba Minch.
#    Site labels are underscore-normalised (e.g. "Arba Minch" -> "Arba_Minch")
#    so they survive shell word-splitting and filename construction.
#    Only monoclonal samples (FWS >= FWS_TH) retained.
# ------------------------------------------------------------------
OUTDIR="tajimasD_by_site"
mkdir -p "$OUTDIR"
echo "=== By Site (Ethiopia 2024) ==="
sites=$(awk -F'\t' -v col_c="$country_col" -v col_y="$year_col" -v col_d="$district_col" \
    -v col_f="$fws_col" -v th="$FWS_TH" \
    'NR>1 && $col_c=="Ethiopia" && $col_y=="2024" && $col_d != "" && $col_f >= th {
        s=$col_d; gsub(/ +/, "_", s); print s
     }' "$METADATA" | sort -u)
for group in $sites; do
    echo "Processing: $group"
    awk -F'\t' -v col_c="$country_col" -v col_y="$year_col" -v col_d="$district_col" \
        -v col_s="$wgs_col" -v col_f="$fws_col" -v th="$FWS_TH" -v g="$group" \
        'NR>1 && $col_c=="Ethiopia" && $col_y=="2024" && $col_f >= th {
            s=$col_d; gsub(/ +/, "_", s)
            if (s == g) print $col_s
         }' "$METADATA" > "${OUTDIR}/${group}_samples.txt"
    run_tajima "$group" "${OUTDIR}/${group}_samples.txt" "$OUTDIR"
done
echo "Sites complete."
echo ""

# ------------------------------------------------------------------
# 8. Sites grouping (Ethiopian 2024 collection sites, collapsed)
#    Co_endemic_sites = Gondar + Mizan + Metehera Fentale, compared with
#    Gambella and Arba Minch.
#    Gambella and Arba_Minch are recomputed here rather than reused from
#    section 7 so this comparison is self-contained; the values are identical
#    to the by_site ones and downstream scripts tell the two apart on the
#    'comparison' column of the merged output.
#    Only monoclonal samples (FWS >= FWS_TH) retained.
# ------------------------------------------------------------------
OUTDIR="tajimasD_sites_grouping"
mkdir -p "$OUTDIR"
echo "=== Sites grouping (Ethiopia 2024) ==="
for group in Co_endemic_sites Gambella Arba_Minch; do
    echo "Processing: $group"
    awk -F'\t' -v col_c="$country_col" -v col_y="$year_col" -v col_d="$district_col" \
        -v col_s="$wgs_col" -v col_f="$fws_col" -v th="$FWS_TH" \
        -v co="$CO_ENDEMIC" -v g="$group" \
        'NR>1 && $col_c=="Ethiopia" && $col_y=="2024" && $col_d != "" && $col_f >= th {
            if ($col_d ~ co) { s = "Co_endemic_sites" }
            else { s = $col_d; gsub(/ +/, "_", s) }
            if (s == g) print $col_s
         }' "$METADATA" > "${OUTDIR}/${group}_samples.txt"
    run_tajima "$group" "${OUTDIR}/${group}_samples.txt" "$OUTDIR"
done
echo "Sites grouping complete."
echo ""

# ------------------------------------------------------------------
# 9. pfcrt K76 vs 76T (co-endemic sites only)
#    Gondar, Mizan and Metehera Fentale, Ethiopia 2024.
#    Groups are defined by K76T (Pf3D7_07_v3:403,625):
#      K76 = 0/0 (wild-type), 76T = 1/1 (mutant)
#    Heterozygous, missing and second-alternate calls are excluded. K76T is
#    multiallelic and therefore absent from the .bi. VCF used elsewhere, so the
#    genotypes come from K76T_longformat.csv; part of that file is phased
#    ("1|1"), so separators are normalised before matching.
#    The groups are named for the codon-76 genotype rather than the
#    CVIET/CVMNK haplotype because M74I and N75E are not variant in this
#    callset — see ehh/pfcrt_genotype_full_pipeline.sh for the detail.
#    K76T_longformat.csv columns: Chrom,Position,Ref,Alt,Sample,Genotype
#    Only monoclonal samples (FWS >= FWS_TH) retained.
# ------------------------------------------------------------------
OUTDIR="tajimasD_pfcrt"
mkdir -p "$OUTDIR"
echo "=== pfcrt: K76 vs 76T (co-endemic sites) ==="
# Alt is quoted and contains a comma ("C,G"), so the CSV is parsed with a
# quote-aware split rather than -F',' as the other lookups do.
python3 -c '
import csv, sys
with open(sys.argv[1], newline="") as fh:
    for row in csv.DictReader(fh):
        gt = row["Genotype"].replace("|", "/")
        if gt == "0/0":
            print(row["Sample"] + "\tK76")
        elif gt == "1/1":
            print(row["Sample"] + "\t76T")
' "$K76T_CSV" > "$OUTDIR/pfcrt_lookup.tmp"
for group in K76 76T; do
    echo "Processing: $group"
    awk -F'\t' -v col_s="$wgs_col" -v col_c="$country_col" -v col_y="$year_col" \
        -v col_d="$district_col" -v col_f="$fws_col" -v th="$FWS_TH" \
        -v co="$CO_ENDEMIC" -v g="$group" \
        'BEGIN{while((getline line < ARGV[1])>0){split(line,a,"\t"); grp[a[1]]=a[2]} ARGV[1]=""}
         NR>1 && $col_c=="Ethiopia" && $col_y=="2024" && $col_d ~ co \
              && ($col_s in grp) && grp[$col_s]==g && $col_f >= th {print $col_s}' \
        "$OUTDIR/pfcrt_lookup.tmp" "$METADATA" > "${OUTDIR}/${group}_samples.txt"
    run_tajima "$group" "${OUTDIR}/${group}_samples.txt" "$OUTDIR"
done
rm "$OUTDIR/pfcrt_lookup.tmp"
echo "pfcrt complete."
echo ""

# ------------------------------------------------------------------
# Step 2: Annotate (bottom 10% windows, gene overlap)
# ------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
echo "=== Annotating ==="
Rscript "$SCRIPT_DIR/annotate.R"
echo ""

# ------------------------------------------------------------------
# Step 3: Merge outputs across all comparisons
# ------------------------------------------------------------------
echo "=== Merging ==="

# 3a. Merge all-windows annotated files — adds a 'comparison' column
ALL_ANNOTATED_MERGED="$DATA_DIR/tajimasD_all_annotated_merged.tsv"
header_written=0
> "$ALL_ANNOTATED_MERGED"
for f in "$DATA_DIR"/tajimasD_annotated_*/*_10kb_annotated.tsv; do
    [[ -f "$f" ]] || continue
    comparison=$(basename "$(dirname "$f")" | sed 's/^tajimasD_annotated_//')
    if [[ $header_written -eq 0 ]]; then
        head -1 "$f" | awk '{print $0 "\tcomparison"}' >> "$ALL_ANNOTATED_MERGED"
        header_written=1
    fi
    tail -n +2 "$f" | awk -v cmp="$comparison" '{print $0 "\t" cmp}' >> "$ALL_ANNOTATED_MERGED"
done
if [[ $header_written -eq 1 ]]; then
    echo "  All-windows annotated merged: $ALL_ANNOTATED_MERGED ($(wc -l < "$ALL_ANNOTATED_MERGED") lines)"
else
    echo "  WARNING: no all-windows annotated files found."
fi

# 3c. Merge annotated q0.1 files — adds a 'comparison' column
ANNOTATED_MERGED="$DATA_DIR/tajimasD_all_annotated_q0.1_merged.tsv"
header_written=0
> "$ANNOTATED_MERGED"
for f in "$DATA_DIR"/tajimasD_annotated_*/*_annotated_q0.1.tsv; do
    [[ -f "$f" ]] || continue
    comparison=$(basename "$(dirname "$f")" | sed 's/^tajimasD_annotated_//')
    if [[ $header_written -eq 0 ]]; then
        head -1 "$f" | awk '{print $0 "\tcomparison"}' >> "$ANNOTATED_MERGED"
        header_written=1
    fi
    tail -n +2 "$f" | awk -v cmp="$comparison" '{print $0 "\t" cmp}' >> "$ANNOTATED_MERGED"
done
if [[ $header_written -eq 1 ]]; then
    echo "  Annotated merged: $ANNOTATED_MERGED ($(wc -l < "$ANNOTATED_MERGED") lines)"
else
    echo "  WARNING: no annotated q0.1 files found."
fi

# 3d. Merge raw Tajima.D files — adds 'group' and 'comparison' columns
RAW_MERGED="$DATA_DIR/tajimasD_all_merged.tsv"
header_written=0
> "$RAW_MERGED"
for dir in tajimasD_by_country tajimasD_by_region tajimasD_old_vs_new \
           tajimasD_north_vs_rest tajimasD_r622i tajimasD_hrp_status \
           tajimasD_by_site tajimasD_sites_grouping tajimasD_pfcrt; do
    [[ -d "$DATA_DIR/$dir" ]] || continue
    comparison=$(echo "$dir" | sed 's/^tajimasD_//')
    for f in "$DATA_DIR/$dir"/*.Tajima.D; do
        [[ -f "$f" ]] || continue
        group=$(basename "$f" | sed 's/^tajimasd_//;s/_10kb\.Tajima\.D$//')
        if [[ $header_written -eq 0 ]]; then
            head -1 "$f" | awk '{print $0 "\tgroup\tcomparison"}' >> "$RAW_MERGED"
            header_written=1
        fi
        tail -n +2 "$f" | awk -v g="$group" -v cmp="$comparison" '{print $0 "\t" g "\t" cmp}' >> "$RAW_MERGED"
    done
done
if [[ $header_written -eq 1 ]]; then
    echo "  Raw merged: $RAW_MERGED ($(wc -l < "$RAW_MERGED") lines)"
else
    echo "  WARNING: no .Tajima.D files found."
fi
echo ""

echo "All done."
